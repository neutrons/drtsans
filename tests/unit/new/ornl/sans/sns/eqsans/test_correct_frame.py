from collections import namedtuple
from os.path import join as pjoin
import pytest
from pytest import approx
from numpy.testing import assert_almost_equal
from mantid.simpleapi import Load, CreateWorkspace
from ornl.sans import wavelength as sans_wavelength
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans import correct_frame
from ornl.settings import amend_config, unique_workspace_name
from ornl.sans.geometry import source_detector_distance

BandsTuple = namedtuple('BandsTuple', 'lead skip')


def test_transmitted_bands(reference_dir):
    with amend_config(data_dir=reference_dir.new.eqsans):
        ws = Load(Filename='EQSANS_86217')
        bands = correct_frame.transmitted_bands(ws)
        assert_almost_equal((bands.lead.min, bands.lead.max),
                            (2.48, 6.78), decimal=2)
        assert_almost_equal((bands.skip.min, bands.skip.max),
                            (10.90, 15.23), decimal=2)


def test_transmitted_bands_clipped(reference_dir):
    with amend_config(data_dir=reference_dir.new.eqsans):
        ws = Load(Filename='EQSANS_86217')
        sdd = source_detector_distance(ws, unit='m')
        bands_0 = correct_frame.transmitted_bands_clipped(ws, sdd, 0.0, 0.0)
        lwc, hwc = (0.139, 0.560)  # expected clippings
        # With no interior clipping
        bands = correct_frame.transmitted_bands_clipped(ws, sdd, 500, 2000,
                                                        interior_clip=False)
        # Check clippings for the lead pulse
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0), 0.01)
        # Check clippings for the skip pulse
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0, b2_0 - hwc), 0.01)
        # With interior clipping
        bands = correct_frame.transmitted_bands_clipped(ws, sdd, 500, 2000,
                                                        interior_clip=True)
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)


@pytest.mark.offline
def test_log_tof_structure(reference_dir):
    file_name = pjoin(reference_dir.new.eqsans, 'test_correct_frame',
                      'EQSANS_92353_no_events.nxs')
    for ny, refv in ((False, 30833), (True, 28333)):
        ws = Load(file_name, OutputWorkspace=unique_workspace_name())
        correct_frame.log_tof_structure(ws, 500, 2000, interior_clip=ny)
        sl = SampleLogs(ws)
        assert sl.tof_frame_width.value == approx(33333, abs=1.0)
        assert sl.tof_frame_width_clipped.value == approx(refv, abs=1.0)
        ws.delete()


def test_band_structure_logs():
    w = CreateWorkspace([0], [0], OutputWorkspace=unique_workspace_name())
    with pytest.raises(RuntimeError, match='Band structure not found in the logs'):
        correct_frame.metadata_bands(w)
    SampleLogs(w).insert('is_frame_skipping', 0)
    correct_frame.log_band_structure(w, BandsTuple(sans_wavelength.Wband(1.5, 2.42), None))
    bands = correct_frame.metadata_bands(w)
    assert bands.lead.min, bands.lead.max == approx(1.5, 2.42)
    SampleLogs(w).insert('is_frame_skipping', 1)
    with pytest.raises(RuntimeError, match='Bands from the skipped pulse missing in the logs'):
        correct_frame.metadata_bands(w)
    correct_frame.log_band_structure(w, BandsTuple(sans_wavelength.Wband(1.5, 2.42),
                                                   sans_wavelength.Wband(6.07, 10.01)))
    bands = correct_frame.metadata_bands(w)
    assert bands.lead.min, bands.lead.max == approx(1.5, 2.42)
    assert bands.skip.min, bands.skip.max == approx(6.07, 10.01)


if __name__ == '__main__':
    pytest.main()
