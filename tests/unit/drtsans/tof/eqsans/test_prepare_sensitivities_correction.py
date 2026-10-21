import json

import numpy as np
import pytest
from mantid.simpleapi import CreateSampleWorkspace, CreateWorkspace, MaskDetectors

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


def test_eqsans_run_inputs_are_normalized_to_one_run():
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs(1)
    preparer.set_direct_beam_runs(["EQSANS_2"])
    preparer.set_dark_current_runs("3")
    preparer.set_blocked_beam_runs(4)
    preparer.set_empty_cell_runs("EQSANS_5.nxs.h5")
    preparer.set_transmission_correction(6, 7, 8)

    assert preparer._flood_runs == [1]
    assert preparer._direct_beam_center_runs == ["EQSANS_2"]
    assert preparer._dark_current_runs == ["3"]
    assert preparer._blocked_beam_runs == [4]
    assert preparer._empty_cell_runs == ["EQSANS_5.nxs.h5"]
    assert preparer._transmission_flood_runs == [6]
    assert preparer._transmission_empty_cell_runs == [7]
    assert preparer._transmission_reference_runs == [8]


def test_eqsans_optional_run_inputs_accept_none():
    preparer = PrepareSensitivityCorrection()
    preparer.set_dark_current_runs(None)
    preparer.set_blocked_beam_runs(None)
    preparer.set_empty_cell_runs(None)
    preparer.set_transmission_correction(None, None, None)

    assert preparer._dark_current_runs is None
    assert preparer._blocked_beam_runs is None
    assert preparer._empty_cell_runs is None
    assert preparer._transmission_reference_runs is None


@pytest.mark.parametrize(
    "setter, name",
    [
        ("set_flood_runs", "flood_runs"),
        ("set_direct_beam_runs", "direct_beam_runs"),
        ("set_dark_current_runs", "dark_current_runs"),
        ("set_blocked_beam_runs", "blocked_beam_runs"),
        ("set_empty_cell_runs", "empty_cell_runs"),
    ],
)
def test_eqsans_run_inputs_reject_more_than_one_run(setter, name):
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match=f"{name} must contain exactly one run"):
        getattr(preparer, setter)([1, 2])


def test_eqsans_transmission_runs_reject_more_than_one_run():
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match="transmission_reference_runs must contain exactly one run"):
        preparer.set_transmission_correction(1, 2, [3, 4])


def test_eqsans_execute_rejects_moving_detector_method(monkeypatch):
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs(1)
    preparer.set_direct_beam_runs(2)
    monkeypatch.setattr(PrepareBase, "execute", lambda *_args, **_kwargs: pytest.fail("common workflow must not run"))

    with pytest.raises(ValueError, match="moving-detector method is not supported for EQSANS"):
        preparer.execute(True, 0.5, 2.0, "sensitivity.nxs")


def test_eqsans_execute_delegates_patch_method_to_common_workflow(monkeypatch):
    preparer = PrepareSensitivityCorrection()
    calls = []
    monkeypatch.setattr(PrepareBase, "execute", lambda *args, **kwargs: calls.append((args[1:], kwargs)))

    preparer.execute(False, 0.5, 2.0, "sensitivity.nxs")

    assert calls == [((False, 0.5, 2.0, "sensitivity.nxs"), {"enforce_use_nexus_idf": False, "debug_mode": False})]


def test_eqsans_transmission_runs_must_be_supplied_together():
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match="must be supplied together"):
        preparer.set_transmission_correction(1, None, 3)


def test_eqsans_tof_cuts_are_independent_clipping_widths():
    preparer = PrepareSensitivityCorrection()

    preparer.set_tof_cuts(low_tof_clip=2000.0, high_tof_clip=500.0)

    assert preparer._low_tof_clip == 2000.0
    assert preparer._high_tof_clip == 500.0


@pytest.mark.parametrize(
    "run, expected",
    [(1236, "EQSANS_1236"), ("1236", "EQSANS_1236"), ("EQSANS_1236", "EQSANS_1236")],
)
def test_eqsans_beam_center_and_parent_runs_use_eqsans_run_names(run, expected):
    preparer = PrepareSensitivityCorrection()

    assert preparer._resolve_beam_center_run(run) == expected
    assert preparer._resolve_parent_flood_run(run) == expected


def test_eqsans_beam_center_and_parent_runs_reject_unknown_names():
    preparer = PrepareSensitivityCorrection()

    with pytest.raises(ValueError, match="direct_beam_run"):
        preparer._resolve_beam_center_run("/missing/run.nxs.h5")
    with pytest.raises(ValueError, match="flood_run"):
        preparer._resolve_parent_flood_run("/missing/run.nxs.h5")


def test_eqsans_execute_accepts_instrument_qualified_run_names(monkeypatch, clean_workspace):
    preparer = PrepareSensitivityCorrection()
    preparer.set_flood_runs("EQSANS_1234")
    preparer.set_direct_beam_runs("EQSANS_1236")

    def create_workspace(name):
        clean_workspace(name)
        clean_workspace(f"{name}_new")
        return CreateWorkspace(DataX=[1.0, 2.0], DataY=[1.0, 1.0], DataE=[0.1, 0.1], NSpec=2, OutputWorkspace=name)

    loaded_runs = []
    beam_center_workspaces = []
    parent_runs = []
    saved_files = []

    def fake_prepare_data(data, output_workspace=None, **_kwargs):
        loaded_runs.append(data)
        return create_workspace(output_workspace or "eqsans_named_run_flood")

    def fake_find_beam_center(workspace):
        beam_center_workspaces.append(str(workspace))
        return 0.01, 0.02, None

    def fake_load_events(run, **_kwargs):
        parent_runs.append(run)
        return create_workspace("eqsans_named_run_parent")

    monkeypatch.setitem(common_sensitivity.PREPARE_DATA, "EQSANS", fake_prepare_data)
    monkeypatch.setitem(common_sensitivity.FIND_BEAM_CENTER, "EQSANS", fake_find_beam_center)
    monkeypatch.setitem(
        common_sensitivity.CALCULATE_SENSITIVITY_CORRECTION,
        common_sensitivity.PATCHING_DETECTORS,
        lambda workspace, **_kwargs: workspace,
    )
    monkeypatch.setattr(eqsans_sensitivity, "NormaliseByCurrent", lambda **kwargs: kwargs["InputWorkspace"])
    monkeypatch.setattr(common_sensitivity, "load_events", fake_load_events)
    monkeypatch.setattr(common_sensitivity, "MaskDetectors", lambda **_kwargs: None)
    monkeypatch.setattr(
        common_sensitivity, "SaveNexusProcessed", lambda **kwargs: saved_files.append(kwargs["Filename"])
    )
    monkeypatch.setattr(preparer, "_get_masked_detectors", lambda _workspace: np.zeros((2, 1), dtype=bool))
    monkeypatch.setattr(preparer, "_mask_beam_center", lambda workspace, _beam_center: workspace)
    monkeypatch.setattr(preparer, "_set_mask_value", lambda workspace, *_args: workspace)

    preparer.execute(False, 0.5, 2.0, "sensitivity.nxs")

    assert loaded_runs == ["EQSANS_1236", "EQSANS_1234"]
    assert beam_center_workspaces == ["BC_EQSANS_1236"]
    assert parent_runs == ["EQSANS_1234"]
    assert saved_files == ["sensitivity.nxs"]


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


def test_eqsans_workspace_preparation_converts_event_workspace_to_histogram(monkeypatch, clean_workspace):
    preparer = PrepareSensitivityCorrection()
    events = CreateSampleWorkspace(
        WorkspaceType="Event",
        NumBanks=1,
        BankPixelWidth=2,
        XUnit="Wavelength",
        XMin=2.0,
        XMax=6.0,
        BinWidth=0.5,
        OutputWorkspace="eqsans_event_flood_test",
    )
    clean_workspace(events)
    # The conversion replaces the workspace in place, invalidating the "events" handle
    expected_spectra = events.getNumberHistograms()
    expected_x = events.readX(0).copy()
    expected_y = events.readY(0).copy()
    expected_e = events.readE(0).copy()

    monkeypatch.setitem(eqsans_sensitivity.PREPARE_DATA, "EQSANS", lambda **_kwargs: events)
    monkeypatch.setattr(eqsans_sensitivity, "NormaliseByCurrent", lambda **kwargs: kwargs["InputWorkspace"])

    workspace = preparer._prepare_eqsans_workspace(123, [0.0, 0.0], None, "flood_run")

    assert workspace.id() == "Workspace2D"
    assert workspace.getNumberHistograms() == expected_spectra
    np.testing.assert_allclose(workspace.readX(0), expected_x)
    np.testing.assert_allclose(workspace.readY(0), expected_y)
    np.testing.assert_allclose(workspace.readE(0), expected_e)


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
    preparer.set_flood_runs(1)
    preparer.set_direct_beam_runs(3)
    preparer.set_dark_current_runs(5)
    preparer.set_blocked_beam_runs(7)
    preparer.set_empty_cell_runs(9)
    preparer.set_transmission_correction(11, 13, 15)
    preparer.set_sample_detector_offsets(12.0, 34.0)
    preparer.set_tof_cuts(500.0, 2000.0)
    preparer.set_beam_center_radius(6.0)
    preparer.set_masks("default.xml", "1-4")
    preparer.scale_components = {"detector1": [1.0, 2.0, 1.0]}
    preparer._last_beam_centers = [(0.1, 0.2)]
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
    assert json.loads(run.getProperty("flood_runs").value) == [1]
    assert run.getProperty("beam_center_x").value == pytest.approx(0.1)
    assert run.getProperty("sample_offset").value == 12.0
    assert json.loads(run.getProperty("scale_components").value) == {"detector1": [1.0, 2.0, 1.0]}
    assert run.getProperty("theta_dependent_transmission").value
    assert run.getProperty("wavelength_weighting").value == "counts"
