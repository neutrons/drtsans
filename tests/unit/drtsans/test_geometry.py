# local imports
from drtsans.instruments import empty_instrument_workspace
from drtsans.mask_utils import apply_mask
from drtsans.settings import unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans import geometry as geo
from drtsans.mono.biosans import geometry as biogeo

# third party imports
from mantid.simpleapi import (
    AddSampleLogMultiple,
    LoadEmptyInstrument,
    LoadEventNexus,
    MoveInstrumentComponent,
)
import numpy as np
from numpy.testing import assert_almost_equal
import pytest

# standard imports
from os.path import join as path_join


@pytest.fixture(scope="module")
def wss():
    r"""Just one workspace for each instrument"""

    # Load an EQSANS instrument and mess with the instrument components
    _eq_ws = LoadEmptyInstrument(InstrumentName="EQSANS")
    geo.sample_detector_distance(_eq_ws)
    for component, shift in (
        ("detector1", 1.3),
        ("sample-position", 0.02),
        ("moderator", 0.3),
    ):
        MoveInstrumentComponent(_eq_ws, ComponentName=component, Z=shift)

    # ssd: source-sample-distance, sdd: sample-detector-distance
    return dict(biosans=None, eqsans=dict(ws=_eq_ws, ssd=13842, sdd=1280), gpsans=None)


def test_source_sample_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.source_sample_distance(v["ws"]) == pytest.approx(v["ssd"], rel=0.01)


def test_sample_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.sample_detector_distance(v["ws"]) == pytest.approx(v["sdd"], rel=0.01)


def test_source_detector_distance(wss):
    for v in wss.values():
        if v is not None:
            assert geo.source_detector_distance(v["ws"]) == pytest.approx(v["ssd"] + v["sdd"], rel=0.01)


def test_detector_translation():
    r"""Ascertain sub-components are moved when main detector is moved"""
    translation = np.array([0.01, 0.1, 1.0])
    detector_name = "detector1"
    for instrument_name in ("EQ-SANS", "CG2"):
        workspace = LoadEmptyInstrument(
            InstrumentName=instrument_name,
            OutputWorkspace=unique_workspace_dundername(),
        )
        instrument = workspace.getInstrument()
        component_detector = instrument.getComponentByName(detector_name)
        component_bank = instrument.getComponentByName("bank42")
        component_detector = instrument.getDetector(42)
        initial_positions = [c.getPos() for c in (component_detector, component_bank, component_detector)]
        MoveInstrumentComponent(
            workspace, ComponentName=detector_name, RelativePosition=True, **dict(zip(("X", "Y", "Z"), translation))
        )
        final_positions = [c.getPos() for c in (component_detector, component_bank, component_detector)]
        for i, final_position in enumerate(final_positions):
            assert final_position == pytest.approx(np.array(initial_positions[i]) + translation, abs=1e-4)
        workspace.delete()


@pytest.mark.parametrize(
    "filename, component, detmin, detmax",
    [
        ("EQ-SANS_Definition.xml", "", 0, 49151),
        ("BIOSANS_Definition.xml", "", 0, 52 * 8 * 256 - 1),
        ("BIOSANS_Definition.xml", "detector1", 0, 24 * 8 * 256 - 1),
        ("BIOSANS_Definition.xml", "wing_detector", 24 * 8 * 256, 44 * 8 * 256 - 1),
        ("BIOSANS_Definition.xml", "midrange_detector", 44 * 8 * 256, 52 * 8 * 256 - 1),
    ],
)
def test_bank_detector_ids(filename, component, detmin, detmax, fetch_idf):
    wksp = LoadEmptyInstrument(Filename=fetch_idf(filename), OutputWorkspace=unique_workspace_dundername())
    num_detectors = detmax - detmin + 1

    # None test
    detIDs = geo.bank_detector_ids(wksp, component=component, masked=None)
    assert detIDs.size == num_detectors
    assert detIDs.min() == detmin
    assert detIDs.max() == detmax

    detIDs = geo.bank_detector_ids(wksp, component=component, masked=False)
    assert detIDs.size == num_detectors
    assert detIDs.min() == detmin
    assert detIDs.max() == detmax

    detIDs = geo.bank_detector_ids(wksp, component=component, masked=True)
    assert len(detIDs) == 0


@pytest.mark.parametrize(
    "filename, component, wksp_index_min, wksp_index_max",
    [
        ("EQ-SANS_Definition.xml", "", 1, 49151 + 2),
        ("BIOSANS_Definition.xml", "", 2, 52 * 8 * 256 + 2),
        ("BIOSANS_Definition.xml", "detector1", 2, 24 * 8 * 256 + 2),
        ("BIOSANS_Definition.xml", "wing_detector", 24 * 8 * 256 + 2, 44 * 8 * 256 + 2),
        ("BIOSANS_Definition.xml", "midrange_detector", 44 * 8 * 256 + 2, 52 * 8 * 256 + 2),
    ],
)
def test_bank_workspace_indices(filename, component, wksp_index_min, wksp_index_max, fetch_idf):
    wksp = LoadEmptyInstrument(Filename=fetch_idf(filename), OutputWorkspace=unique_workspace_dundername())

    wksp_indices = geo.bank_workspace_index_range(wksp, component)
    assert wksp_indices[0] >= 0
    assert wksp_indices[1] <= wksp.getNumberHistograms()
    assert wksp_indices[0] == wksp_index_min
    assert wksp_indices[1] == wksp_index_max


@pytest.mark.parametrize(
    "workspace_with_instrument",
    [
        {
            "instrument_geometry": "n-pack",
            "n_tubes": 2,
            "n_pixels": 2,
            "spacing": 0.0,
            "x_center": 0.0,
            "y_center": 0.0,
            "z_center": 0.0,  # detector center
            "diameter": 0.02,
            "height": 0.02,  # pixel dimensions
            # ruff keeps removing the following 3 lines
            # "x_center": 0.0,
            # "y_center": 0.0,
            # "z_center": 0.0,
        }
    ],
    indirect=True,
)
def test_pixel_centers(workspace_with_instrument):
    r"""
    Pixel centers for a detector array made up of two tubes, each with two pixels.
    There's no spacing between tubes.
    Detector is 1 meter away from the sample
    The shape of a detector pixel is a cylinder of 20mm diameter and 20mm in height.
    """
    # The generated IDF file still puts pixel position at the bottome edge center
    input_workspace = unique_workspace_dundername()  # temporary workspace
    workspace_with_instrument(
        axis_values=[2.0, 2.1],
        intensities=[[1.0, 2.0], [3.0, 2.0]],
        output_workspace=input_workspace,
    )
    pixel_positions = geo.pixel_centers(input_workspace, [0, 1, 2, 3])
    expected = 1.0e-03 * np.array([[10, -20, 0.0], [10, 0, 0.0], [-10, -20, 0.0], [-10, 0, 0.0]])  # in meters
    assert pixel_positions == pytest.approx(expected)


def test_logged_smearing_pixel_size(workspace_with_instrument):
    workspace = workspace_with_instrument()

    # Assert default value of `None` when no smearing pixels are provided
    logged_values = geo.logged_smearing_pixel_size(workspace)
    assert list(logged_values) == [None, None]

    # Assert values are correctly retrieved when smearing pixels are provided
    values, names, units = (
        [0.0042, 0.0024],
        ["smearingPixelSizeX", "smearingPixelSizeY"],
        ["m", "m"],
    )
    AddSampleLogMultiple(Workspace=workspace, LogNames=names, LogValues=values, LogUnits=units)
    # the order of `logged_values` should be the same as that of `values`
    logged_values = geo.logged_smearing_pixel_size(workspace)
    assert logged_values == pytest.approx(logged_values)


def test_sample_aperture_diameter(serve_events_workspace, reference_dir):
    # NOTE:
    # serve_events_workspace is a pytest fixture that is hardcoded to work for
    # eqsans tests.
    input_workspace = serve_events_workspace("EQSANS_92353.nxs.h5")
    # diameter is retrieved from log 'beamslit4', and we convert the 10mm into 0.01 meters
    assert geo.sample_aperture_diameter(input_workspace, unit="m") == pytest.approx(0.01, abs=0.1)
    # verify entry 'sample_aperture_diameter' has been added to the logs
    assert SampleLogs(input_workspace).single_value("sample_aperture_diameter") == pytest.approx(10.0, abs=0.1)
    # test a run containing "sample_aperture_radius" instead of "sample_aperture_diameter"
    workspace = LoadEventNexus(
        Filename=path_join(reference_dir.gpsans, "geometry", "CG2_1338.nxs.h5"),
        OutputWorkspace=unique_workspace_dundername(),
        MetaDataOnly=True,
        LoadLogs=True,
    )
    assert geo.sample_aperture_diameter(workspace, unit="mm") == pytest.approx(14.0, abs=0.1)
    workspace.delete()


def test_source_aperture_diameter(reference_dir):
    # test a run containing "sample_aperture_radius" instead of "sample_aperture_diameter"
    workspace = LoadEventNexus(
        Filename=path_join(reference_dir.gpsans, "geometry", "CG2_1338.nxs.h5"),
        OutputWorkspace=unique_workspace_dundername(),
        MetaDataOnly=True,
        LoadLogs=True,
    )
    assert geo.source_aperture_diameter(workspace, unit="mm") == pytest.approx(40.0, abs=0.1)
    workspace.delete()


def test_translate_source_by_z(reference_dir):
    filename = path_join(reference_dir.gpsans, "geometry", "CG2_1338.nxs.h5")
    workspace = LoadEventNexus(
        Filename=filename,
        OutputWorkspace=unique_workspace_dundername(),
        MetaDataOnly=False,
        LoadLogs=True,
    )
    geo.translate_source_by_z(workspace)
    assert workspace.getInstrument().getComponentByName("moderator").getPos().Z() == pytest.approx(-7.283, abs=0.1)


@pytest.mark.parametrize(
    "instrument, pixel_size",
    [
        ("BIOSANS", (0.00804, 0.00409)),
        ("EQSANS", (0.00804, 0.00409)),
        ("GPSANS", (0.00804, 0.00409)),
    ],
)
def test_nominal_pixel_size(instrument, pixel_size):
    workspace = LoadEmptyInstrument(InstrumentName=instrument, OutputWorkspace=unique_workspace_dundername())
    assert geo.nominal_pixel_size(workspace) == pytest.approx(pixel_size, abs=1.0e-04)
    workspace.delete()


def test_get_position_south_detector(fetch_idf, temp_workspace_name):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename=fetch_idf("BIOSANS_Definition.xml"), OutputWorkspace=temp_workspace_name()
    )
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=7.00, RelativePosition=False)
    assert_almost_equal(geo.get_position_south_detector(workspace), 7.000, decimal=3)


def test_get_curvature_radius(temp_workspace_name):
    workspace = LoadEmptyInstrument(
        InstrumentName="BIOSANS", Filename="BIOSANS_Definition.xml", OutputWorkspace=temp_workspace_name()
    )
    assert_almost_equal(geo.get_curvature_radius(workspace, "wing_detector"), 1.1633, decimal=4)
    assert_almost_equal(geo.get_curvature_radius(workspace, "midrange_detector"), 4.0000, decimal=4)


def test_pixel_masks(temp_workspace_name):
    workspace = empty_instrument_workspace(temp_workspace_name(), filename="BIOSANS_Definition.xml")
    apply_mask(workspace, Components="wing_detector")
    assert np.alltrue(geo.get_pixel_masks(workspace, "wing_detector"))
    apply_mask(workspace, mask=list(range(10)))
    mask_in_detector1 = geo.get_pixel_masks(workspace, "wing_detector")[:10]
    assert np.alltrue(mask_in_detector1)


def test_get_pixel_distances(temp_workspace_name, fetch_idf):
    workspace = empty_instrument_workspace(temp_workspace_name(), filename=fetch_idf("BIOSANS_Definition.xml"))
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=7.00, RelativePosition=False)
    expected_for_component = {
        "detector1": (7.0358, 7.0439),
        "wing_detector": (1.2402, 1.2477),
        "midrange_detector": (3.9962, 4.0043),
    }
    for component, expected in expected_for_component.items():
        distances = geo.get_pixel_distances(workspace, component)
        assert_almost_equal((distances[0], distances[-1]), expected, decimal=3)


def test_get_twothetas(temp_workspace_name, fetch_idf):
    workspace = empty_instrument_workspace(temp_workspace_name(), filename=fetch_idf("BIOSANS_Definition.xml"))
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=7.00, RelativePosition=False)
    expected_for_component = {
        "detector1": (6.105, 6.098),
        "wing_detector": (24.836, 49.498),
        "midrange_detector": (8.980, 7.475),
    }
    for component, expected in expected_for_component.items():
        twothetas = geo.get_twothetas(workspace, component, units="degrees")
        assert_almost_equal((twothetas[0], twothetas[-1]), expected, decimal=2)
        x, y, z = geo.get_xyz(workspace, component)
        twothetas_other = np.degrees(np.arctan(np.sqrt(x**2 + y**2) / z))
        assert_almost_equal(twothetas, twothetas_other, decimal=1)


def test_get_solid_angles(temp_workspace_name, fetch_idf):
    workspace = empty_instrument_workspace(temp_workspace_name(), filename=fetch_idf("BIOSANS_Definition.xml"))
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=7.00, RelativePosition=False)
    expected_for_component = {
        "detector1": (0.02009, 0.0100),
        "wing_detector": (0.6097, 0.2994),
        "midrange_detector": (0.0626, 0.0311),
    }
    for component, expected in expected_for_component.items():
        solid_angles = geo.get_solid_angles(workspace, component)
        assert_almost_equal((solid_angles[0], solid_angles[-1]), expected, decimal=3)
    #
    # move all components. Only solid angles for detector1 should change
    MoveInstrumentComponent(Workspace=workspace, ComponentName="detector1", Z=5.00, RelativePosition=False)
    biogeo.set_angle_wing_detector(workspace, angle=30.0)  # degrees
    biogeo.set_angle_midrange_detector(workspace, angle=5.0)  # degrees
    expected_for_component["detector1"] = (0.03878, 0.0192)
    for component, expected in expected_for_component.items():
        solid_angles = geo.get_solid_angles(workspace, component, back_panel_attenuation=0.5)
        assert_almost_equal((solid_angles[0], solid_angles[-1]), expected, decimal=3)


if __name__ == "__main__":
    pytest.main([__file__])
