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

# public API
from ornl.sans.sns import eqsans
from ornl.sans import solid_angle_correction

# protected API
from ornl.settings import (amend_config, namedtuplefy,
                           unique_workspace_dundername as uwd)
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.geometry import detector_name


keys = ('run', 'num_events', 'sdd', 'ssd', 'min_tof', 'max_tof',
        'skip_frame', 'w_min', 'w_max', 'flux')

values = (('EQSANS_86217', 508339, 1300, 14122, 9717, 59764, True,
           2.61, 14.72, 1654),
          ('EQSANS_92353', 262291, 4000, 14122, 11410, 61439, True,
           2.59, 12.98, 1299),
          ('EQSANS_85550', 270022, 5000, 14122, 12036, 62060, True,
           2.59, 12.43, 1804),
          ('EQSANS_101595', 289989, 1300, 14122, 7781, 24511, False,
           2.11, 5.65, 429),
          ('EQSANS_88565', 19362, 4000, 14122, 45615, 62302, False,
           10.02, 13.2, 2143),
          ('EQSANS_88901', 340431, 8000, 14122, 67332, 83998, False,
           11.99, 14.62, 157523))

run_sets = [{k: v for k, v in zip(keys, value)} for value in values]


@pytest.fixture(scope='module')
def flux_file(refd):
    return pj(refd.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt')


@pytest.fixture(scope='module', params=run_sets)
@namedtuplefy
def rs(refd, request):
    run_set = request.param
    run = run_set['run']
    ws = uwd()
    with amend_config(data_dir=refd.new.eqsans):
        eqsans.load_events(run, output_workspace=ws)
    kw = dict(low_tof_clip=500, high_tof_clip=2000, output_workspace=uwd())
    wl = eqsans.transform_to_wavelength(ws, **kw)
    return {**run_set, **dict(ws=ws, wl=wl)}


class TestLoadEvents(object):

    def test_loading_file(self, rs):
        ws = mtd[rs.ws]
        assert isinstance(ws, EventWorkspace)
        assert ws.getNumberEvents() == rs.num_events

    def test_geometry(self, rs):
        ws = mtd[rs.ws]
        # assert distance of detector1 same as that in detectorZ of the logs
        instrument = ws.getInstrument()
        det = instrument.getComponentByName(detector_name(instrument))
        d1 = det.getDistance(instrument.getSample())
        assert rs.sdd == pytest.approx(d1 * 1000, abs=1)
        # Check logs
        sl = SampleLogs(ws)
        assert rs.ssd == approx(sl.single_value('source-sample-distance'),
                                abs=1)
        assert rs.sdd == approx(sl.single_value('sample-detector-distance'),
                                abs=1)

    def test_tofs(self, rs):
        ws = mtd[rs.ws]
        assert ws.getTofMin() == pytest.approx(rs.min_tof, abs=1)
        assert ws.getTofMax() == pytest.approx(rs.max_tof, abs=1)
        assert bool(SampleLogs(ws).is_frame_skipping.value) == rs.skip_frame

    def test_offsets(self, refd):
        with amend_config(data_dir=refd.new.eqsans):
            ws = eqsans.load_events('EQSANS_86217', output_workspace=uwd(),
                                    detector_offset=42, sample_offset=24)
        sl = SampleLogs(ws)
        ssd = sl.single_value('source-sample-distance')
        sdd = sl.single_value('sample-detector-distance')
        assert 14146 == pytest.approx(ssd, abs=0.1)
        assert 1318 == pytest.approx(sdd, abs=0.1)


def test_transform_to_wavelength(rs):
    ws = eqsans.transform_to_wavelength(rs.ws, low_tof_clip=500,
                                        high_tof_clip=2000,
                                        zero_uncertainty=42.0,
                                        output_workspace=uwd())
    sl = SampleLogs(ws)
    assert sl.wavelength_min.value == approx(rs.w_min, abs=0.05)
    assert sl.wavelength_max.value == approx(rs.w_max, abs=0.05)
    # assert zero uncertainty assignment
    for i in range(ws.getNumberHistograms()):
        zci = np.where(ws.dataY(i) == 0)[0]  # zero count indices
        assert ws.dataE(i)[zci] == approx(42.0*np.ones(len(zci)))


def test_normalise_by_flux(rs, flux_file):
    ws = eqsans.normalise_by_flux(rs.wl, flux_file, output_workspace=uwd())
    ws = SumSpectra(ws)
    assert np.average(ws.dataY(0)) == approx(rs.flux, abs=1)


def test_subtract_background(refd):
    data_dir = pj(refd.new.eqsans, 'test_subtract_background')
    ws = LoadNexus(pj(data_dir, 'sample.nxs'), OutputWorkspace=uwd())
    ws_name = ws.name()
    wb = LoadNexus(pj(data_dir, 'background.nxs'), OutputWorkspace=uwd())
    ws_wb = eqsans.subtract_background(ws, wb, scale=0.42)
    assert ws_wb.name() == ws_name
    assert max(ws_wb.dataY(0)) < 1.e-09


def test_prepare_monitors(refd):
    data_dirs = [refd.new.eqsans, pj(refd.new.eqsans, 'test_integration_api')]
    with amend_config(data_dir=data_dirs):
        # Raises for a run in skip frame mode
        with pytest.raises(RuntimeError, match='cannot correct monitor'):
            eqsans.prepare_monitors('EQSANS_92353')
        # Compare for a run in non-skip frame mode
        w = eqsans.prepare_monitors('EQSANS_88565')
        assert w.name() == 'EQSANS_88565_monitors'
        sl = SampleLogs(w)
        assert sl.wavelength_min.value == approx(9.9, abs=0.1)
        assert sl.wavelength_max.value == approx(13.6, abs=0.1)
        v = LoadNexus('EQSANS_88565_monitors_wav.nxs', OutputWorkspace=uwd())
        assert CompareWorkspaces(w, v, Tolerance=1e-3,
                                 ToleranceRelErr=True).Result is True


@pytest.mark.skip(reason="prepare data not yet completed")
def test_prepared_data(eqsans_f):
    """
        This is Section 3 of the requirements document.
        It should just be a convenience function that brings together
        loading, moving the detector, normalizing, binning in wavelength,
        and subtracting dark current.
    """
    ws = eqsans.prepare_data(eqsans_f['data'])
    assert isinstance(ws, EventWorkspace)


def test_solid_angle(rs):
    ws2 = solid_angle_correction(rs.ws, output_workspace=uwd(),
                                 detector_type='VerticalTube')
    assert isinstance(ws2, EventWorkspace)
    assert ws2.getNumberEvents() == rs.num_events


@pytest.mark.parametrize('name', ['save_ascii_1D', 'save_xml_1D',
                                  'save_nist_dat', 'save_nexus'])
def test_api_contents(name):
    assert name in dir(eqsans)


if __name__ == '__main__':
    pytest.main()
