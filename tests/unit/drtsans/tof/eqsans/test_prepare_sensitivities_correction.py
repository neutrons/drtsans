import json

import numpy as np
import pytest
from mantid.simpleapi import CreateWorkspace, MaskDetectors

from drtsans import prepare_sensivities_correction as common_sensitivity
from drtsans.prepare_sensivities_correction import PrepareSensitivityCorrection as PrepareBase
from drtsans.tof.eqsans import prepare_sensitivities_correction as eqsans_sensitivity
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


def test_eqsans_workspace_preparation_uses_charge_only_normalization(monkeypatch):
    preparer = PrepareSensitivityCorrection()
    preparer.set_sample_detector_offsets(sample_offset=12.0, detector_offset=34.0)
    preparer.set_tof_cuts(low_tof_clip=1650.0, high_tof_clip=3150.0)

    class FakeWorkspace:
        def __str__(self):
            return "flood_workspace"

    prepared = []
    normalized = []

    def fake_prepare_data(**kwargs):
        prepared.append(kwargs)
        return FakeWorkspace()

    def fake_normalise_by_current(**kwargs):
        normalized.append(kwargs)
        return kwargs["InputWorkspace"]

    monkeypatch.setitem(eqsans_sensitivity.PREPARE_DATA, "EQSANS", fake_prepare_data)
    monkeypatch.setattr(eqsans_sensitivity, "NormaliseByCurrent", fake_normalise_by_current)

    workspace = preparer._prepare_eqsans_workspace(123, [0.1, 0.2], None, "flood_run")

    assert str(workspace) == "flood_workspace"
    assert prepared[0]["data"] == "EQSANS_123"
    assert prepared[0]["flux_method"] is None
    assert prepared[0]["sample_offset"] == 12.0
    assert prepared[0]["detector_offset"] == 34.0
    assert prepared[0]["low_tof_clip"] == 1650.0
    assert prepared[0]["high_tof_clip"] == 3150.0
    assert normalized[0]["InputWorkspace"] is workspace
    assert normalized[0]["OutputWorkspace"] == "flood_workspace"


def test_eqsans_flood_subtracts_blocked_beam_and_empty_cell_before_integration(monkeypatch):
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs(1)
    preparer.set_blocked_beam_runs(2)
    preparer.set_empty_cell_runs(3)

    class FakeWorkspace:
        def __init__(self, name):
            self.name = name

        def __str__(self):
            return self.name

        def blocksize(self):
            return 2

    prepared = []
    subtractions = []
    integrations = []

    def fake_prepare_workspace(run, beam_center, dark_current_run, name):
        prepared.append((run, beam_center, dark_current_run, name))
        return FakeWorkspace(name)

    def fake_subtract_background(input_workspace, background, output_workspace=None):
        subtractions.append((input_workspace.name, background.name, output_workspace))
        return input_workspace

    def fake_integration(**kwargs):
        integrations.append(kwargs)
        return kwargs["InputWorkspace"]

    monkeypatch.setattr(preparer, "_prepare_eqsans_workspace", fake_prepare_workspace)
    monkeypatch.setattr(eqsans_sensitivity, "subtract_background", fake_subtract_background)
    monkeypatch.setattr(eqsans_sensitivity, "Integration", fake_integration)

    result = preparer._prepare_flood_data(1, [0.1, 0.2], "EQSANS_10")

    assert result.name == "flood_run"
    assert [run for run, _, _, _ in prepared] == [1, 2, 3]
    assert all(dark_current_run == "EQSANS_10" for _, _, dark_current_run, _ in prepared)
    assert [(left, right) for left, right, _ in subtractions] == [
        ("flood_run", "blocked_beam_run"),
        ("empty_cell_run", "blocked_beam_run"),
        ("flood_run", "empty_cell_run"),
    ]
    assert len(integrations) == 1


def test_eqsans_workspace_preparation_propagates_missing_proton_charge(monkeypatch):
    preparer = PrepareSensitivityCorrection()

    class FakeWorkspace:
        def __str__(self):
            return "flood_workspace"

    def fake_prepare_data(**kwargs):
        return FakeWorkspace()

    def fake_normalise_by_current(**kwargs):
        raise RuntimeError("proton charge log is missing")

    monkeypatch.setitem(eqsans_sensitivity.PREPARE_DATA, "EQSANS", fake_prepare_data)
    monkeypatch.setattr(eqsans_sensitivity, "NormaliseByCurrent", fake_normalise_by_current)

    with pytest.raises(RuntimeError, match="proton charge"):
        preparer._prepare_eqsans_workspace(123, [0.0, 0.0], None, "flood_run")


def test_eqsans_transmission_is_applied_to_flood_and_empty_cell_before_subtraction(monkeypatch):
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs(1)
    preparer.set_blocked_beam_runs(2)
    preparer.set_empty_cell_runs(3)
    preparer.set_transmission_correction(4, 5, 6)

    class FakeWorkspace:
        def __init__(self, name):
            self.name = name

        def __str__(self):
            return self.name

        def blocksize(self):
            return 2

    prepared = []
    subtractions = []
    transmissions = []
    calculated_transmissions = []
    integrations = []

    def fake_prepare_workspace(run, beam_center, dark_current_run, name, solid_angle=None):
        prepared.append((run, name, solid_angle))
        return FakeWorkspace(name)

    def fake_subtract_background(input_workspace, background, output_workspace=None):
        subtractions.append((input_workspace.name, background.name))
        return input_workspace

    def fake_calculate_transmission(sample, reference):
        calculated_transmissions.append((sample.name, reference.name))
        return FakeWorkspace(f"{sample.name}_transmission")

    def fake_apply_transmission(input_workspace, trans_workspace, theta_dependent, output_workspace):
        transmissions.append((input_workspace.name, trans_workspace.name, theta_dependent, output_workspace))
        return input_workspace

    def fake_integration(**kwargs):
        integrations.append(kwargs)
        return kwargs["InputWorkspace"]

    monkeypatch.setattr(preparer, "_prepare_eqsans_workspace", fake_prepare_workspace)
    monkeypatch.setattr(eqsans_sensitivity, "subtract_background", fake_subtract_background)
    monkeypatch.setattr(eqsans_sensitivity, "calculate_transmission", fake_calculate_transmission)
    monkeypatch.setattr(eqsans_sensitivity, "apply_transmission_correction", fake_apply_transmission)
    monkeypatch.setattr(eqsans_sensitivity, "Integration", fake_integration)

    result = preparer._prepare_flood_data(1, [0.1, 0.2], None)

    assert result.name == "flood_run"
    assert [run for run, _, _ in prepared] == [1, 2, 3, 6, 4, 5]
    assert [solid_angle for _, _, solid_angle in prepared] == [None, None, None, False, False, False]
    assert subtractions == [
        ("flood_run", "blocked_beam_run"),
        ("empty_cell_run", "blocked_beam_run"),
        ("flood_run", "empty_cell_run"),
    ]
    assert calculated_transmissions == [
        ("transmission_flood_run", "transmission_reference_run"),
        ("transmission_empty_cell_run", "transmission_reference_run"),
    ]
    assert transmissions == [
        ("flood_run", "transmission_flood_run_transmission", True, "flood_run"),
        ("empty_cell_run", "transmission_empty_cell_run_transmission", True, "empty_cell_run"),
    ]
    assert len(integrations) == 1


def test_eqsans_transmission_requires_an_empty_cell_run():
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs(1)
    preparer.set_transmission_correction(2, 3, 4)

    with pytest.raises(ValueError, match="empty_cell_runs"):
        preparer._prepare_flood_data(1, [0.0, 0.0], None)


@pytest.mark.parametrize("workspace_with_instrument", [dict(name="EQSANS", Nx=2, Ny=1)], indirect=True)
def test_mask_detection_uses_detector_mask_state_not_uncertainty(workspace_with_instrument, clean_workspace):
    workspace = workspace_with_instrument(
        axis_values=[1.0, 2.0],
        intensities=[[1.0], [1.0]],
        uncertainties=[[1.0e-8], [1.0e-8]],
    )
    clean_workspace(workspace)
    MaskDetectors(Workspace=workspace, WorkspaceIndexList=[1])

    masked = PrepareBase._get_masked_detectors(workspace)

    np.testing.assert_array_equal(masked, [[False], [True]])


@pytest.mark.parametrize("workspace_with_instrument", [dict(name="EQSANS", Nx=3, Ny=1)], indirect=True)
def test_set_mask_value_preserves_unmasked_uncertainty(workspace_with_instrument, clean_workspace):
    workspace = workspace_with_instrument(
        axis_values=[1.0, 2.0],
        intensities=[[1.0], [1.0], [1.0]],
        uncertainties=[[0.25], [0.5], [0.75]],
    )
    clean_workspace(workspace)
    MaskDetectors(Workspace=workspace, WorkspaceIndexList=[0, 1])

    PrepareBase._set_mask_value(workspace, np.array([[True], [False], [False]]), False)

    assert np.isnan(workspace.y(0)[0])
    assert np.isnan(workspace.e(0)[0])
    assert np.isneginf(workspace.y(1)[0])
    assert np.isneginf(workspace.e(1)[0])
    assert workspace.y(2)[0] == 1.0
    assert workspace.e(2)[0] == 0.75


def test_masking_beam_center_does_not_reset_propagated_uncertainties(monkeypatch, clean_workspace):
    workspace = CreateWorkspace(
        DataX=[1.0, 2.0],
        DataY=[2.0, 3.0],
        DataE=[0.2, 0.3],
        NSpec=2,
        OutputWorkspace="eqsans_beam_center_uncertainty_test",
    )
    clean_workspace(workspace)
    preparer = PrepareSensitivityCorrection()
    monkeypatch.setattr(common_sensitivity, "circular_mask_from_beam_center", lambda *_: [])
    monkeypatch.setattr(common_sensitivity, "apply_mask", lambda *_args, **_kwargs: None)

    result = preparer._mask_beam_center(workspace, (0.0, 0.0))

    assert result is workspace
    np.testing.assert_allclose(result.extractE().flatten(), [0.2, 0.3])


def test_eqsans_sensitivity_provenance_logs(clean_workspace):
    workspace = CreateWorkspace(
        DataX=[1.0, 2.0],
        DataY=[1.0],
        DataE=[0.1],
        NSpec=1,
        OutputWorkspace="eqsans_provenance_test",
    )
    flood_workspace = CreateWorkspace(
        DataX=[1.0, 2.0],
        DataY=[1.0],
        DataE=[0.1],
        NSpec=1,
        OutputWorkspace="eqsans_provenance_flood_test",
    )
    clean_workspace(workspace)
    clean_workspace(flood_workspace)
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs([1, 2])
    preparer.set_direct_beam_runs([3, 4])
    preparer.set_dark_current_runs([5, 6])
    preparer.set_blocked_beam_runs([7, 8])
    preparer.set_empty_cell_runs([9, 10])
    preparer.set_transmission_correction([11, 12], [13, 14], [15, 16])
    preparer.set_sample_detector_offsets(12.0, 34.0)
    preparer.set_tof_cuts(500.0, 2000.0)
    preparer.set_beam_center_radius(6.0)
    preparer.set_masks("default.xml", "1-4")
    preparer.scale_components = {"detector1": [1.0, 2.0, 1.0]}
    preparer._last_beam_centers = [(0.1, 0.2), (0.3, 0.4)]
    preparer._last_flood_workspace = flood_workspace
    preparer._last_sensitivity_options = {
        "use_moving_detector_method": False,
        "min_threshold": 0.5,
        "max_threshold": 2.0,
        "polynomial_order": 2,
        "min_detectors_per_tube": 50,
    }

    preparer._add_sensitivity_provenance(workspace)

    run = workspace.getRun()
    assert run.getProperty("sensitivity_method").value == "patch"
    assert json.loads(run.getProperty("flood_runs").value) == [1, 2]
    assert run.getProperty("beam_center_x").value == "[0.1, 0.3]"
    assert run.getProperty("sample_offset").value == 12.0
    assert json.loads(run.getProperty("scale_components").value) == {"detector1": [1.0, 2.0, 1.0]}
    assert run.getProperty("theta_dependent_transmission").value
    assert run.getProperty("wavelength_weighting").value == "counts"
