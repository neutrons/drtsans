import pytest

from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection as PrepareBase
from drtsans.tof.eqsans.prepare_sensitivities_correction import PrepareSensitivityCorrection


def test_eqsans_options_are_not_added_to_common_preparer_interface():
    assert not hasattr(PrepareBase, "set_sample_detector_offsets")
    assert not hasattr(PrepareBase, "set_tof_cuts")


def test_eqsans_prepare_data_options_include_configured_geometry_and_tof():
    preparer = PrepareSensitivityCorrection()
    preparer.set_sample_detector_offsets(sample_offset=314.5, detector_offset=80.0)
    preparer.set_tof_cuts(low_tof_clip=1650.0, high_tof_clip=3150.0)

    options = preparer._prepare_data_opts([0.1, 0.2])

    assert options["center_x"] == 0.1
    assert options["center_y"] == 0.2
    assert options["sample_offset"] == 314.5
    assert options["detector_offset"] == 80.0
    assert options["low_tof_clip"] == 1650.0
    assert options["high_tof_clip"] == 3150.0


def test_eqsans_beam_center_options_omit_flood_geometry_and_tof():
    preparer = PrepareSensitivityCorrection()
    preparer.set_sample_detector_offsets(sample_offset=314.5, detector_offset=80.0)
    preparer.set_tof_cuts(low_tof_clip=1650.0, high_tof_clip=3150.0)

    options = preparer._prepare_beam_center_data_opts([0.0, 0.0])

    assert options == {"center_x": 0.0, "center_y": 0.0, "flux_method": None}


def test_eqsans_auxiliary_run_inputs_are_normalized_and_broadcast():
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs([1, 2])
    preparer.set_blocked_beam_runs(3)
    preparer.set_empty_cell_runs("EQSANS_4.nxs.h5")
    preparer.set_transmission_correction(5, 6, 7)

    assert preparer._blocked_beam_runs == [3, 3]
    assert preparer._empty_cell_runs == ["EQSANS_4.nxs.h5", "EQSANS_4.nxs.h5"]
    assert preparer._transmission_flood_runs == [5, 5]
    assert preparer._transmission_empty_cell_runs == [6, 6]
    assert preparer._transmission_reference_runs == [7, 7]


def test_eqsans_auxiliary_run_lists_must_match_flood_runs():
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs([1, 2])

    with pytest.raises(ValueError, match="blocked_beam_runs"):
        preparer.set_blocked_beam_runs([3, 4, 5])


def test_eqsans_transmission_runs_must_be_supplied_together():
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match="must be supplied together"):
        preparer.set_transmission_correction(1, None, 3)


def test_eqsans_tof_cuts_must_be_non_negative():
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match="non-negative"):
        preparer.set_tof_cuts(low_tof_clip=-1.0, high_tof_clip=2000.0)
