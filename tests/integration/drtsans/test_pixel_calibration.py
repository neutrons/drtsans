# local imports
from drtsans.mono.biosans.simulated_intensities import clone_component_intensities
from drtsans.pixel_calibration import (
    BarPositionFormula,
    as_intensities,
    calculate_apparent_tube_width,
    calculate_barscan_calibration,
    find_edges,
    fit_positions,
    load_calibration,
)
from drtsans.samplelogs import SampleLogs
from drtsans.settings import namedtuplefy
from drtsans.tubecollection import TubeCollection

# third-party imports
from mantid.simpleapi import (
    AddSampleLog,
    DeleteWorkspace,
    DeleteWorkspaces,
    LoadEmptyInstrument,
    LoadEventNexus,
    LoadNexus,
    LoadNexusProcessed,
    SaveNexus,
    mtd,
)
import numpy as np

# standard imports
from os.path import join as path_join
from os import listdir
import pytest
import random
import tempfile


def _linear_density(workspace, component="detector1"):
    r"""Tube total intensity per unit length of tube width"""
    collection = TubeCollection(workspace, component).sorted(view="fbfb")
    intensities = np.array([np.sum(tube.readY) for tube in collection])
    widths = np.array([tube[0].width for tube in collection])
    return list(intensities / widths)


def _amplitude(density):
    r"""ratio of fluctuations to the mean, a sort of amplitude if ``density`` is wave-like"""
    return np.std(density) / np.mean(density)


def test_find_edges():
    r"""Finding the edges of the barscan in a single tube,
    then calculate the position and width of the pixels

    based on BarScanShadow_test_KCL_SVP.xlsx and BarScanFitTest_KCL.xlsx
    Testing appendix 2.1 in the master document

    devs - Andrei Savici <saviciat@ornl.gov>,
           Jose Borreguero <borreguerojm@ornl.gov>
    SME  - Ken Littrell <littrellkc@ornl.gov>
    """
    # tube pixel_intensities
    intensities = np.array(
        [
            2,
            2,
            38,
            38,
            38,
            34,
            38,
            41,
            35,
            3,
            4,
            3,
            3,
            4,
            30,
            30,
            37,
            33,
            31,
            39,
            42,
            42,
            2,
            2,
            2,
        ]
    )
    # find edges
    edges = find_edges(intensities)
    # check expected values
    assert (edges.bottom_pixel, edges.top_pixel) == (2, 21)
    assert (edges.bottom_shadow_pixel, edges.above_shadow_pixel) == (9, 14)
    # check if the algorithm fails if there are not enough illuminated pixels
    with pytest.raises(RuntimeError, match="Faulty tube found"):
        find_edges(intensities, min_illuminated_length=len(intensities))


def test_no_shadow():
    r"""Check the failure mode for not finding shaddow"""
    intensities = [1] * 25  # all pixel_intensities are the same, no shaddow
    with pytest.raises(IndexError, match="Could not find bottom shadow edge"):
        find_edges(intensities)


def test_no_bottom_tube_edge():
    r"""Check the failure mode for not finding tube edge pixels"""
    intensities = np.ones(25)
    intensities[14] = 3000  # all pixels but one are below the threshold
    with pytest.raises(IndexError, match="Could not find bottom tube edge"):
        find_edges(intensities)


def test_fit_positions():
    r"""Test the fitted positions and heights of pixels
    Using drtsans.barscan.fit_positions
    """
    # input edge pixels and positions
    edge_pixel = np.arange(25) + 1
    pos = [
        -120,
        -110,
        -98,
        -87,
        -74,
        -62,
        -49,
        -36,
        -23,
        -11,
        1,
        12,
        23,
        32,
        41,
        49,
        57,
        64,
        71,
        78,
        86,
        93,
        102,
        110,
        120,
    ]
    # fit the positions
    fit_results = fit_positions(edge_pixel, pos, tube_pixels=26)

    # compare to the expected data
    expected_positions = [
        -119.659,
        -109.907,
        -98.892,
        -86.953,
        -74.396,
        -61.496,
        -48.494,
        -35.599,
        -22.989,
        -10.808,
        0.8303,
        11.846,
        22.189,
        31.846,
        40.831,
        49.193,
        57.013,
        64.403,
        71.508,
        78.506,
        85.606,
        93.049,
        101.109,
        110.094,
        120.340,
    ]
    expected_heights = [
        9.001,
        10.443,
        11.530,
        12.296,
        12.771,
        12.989,
        12.981,
        12.779,
        12.417,
        11.926,
        11.337,
        10.685,
        10.000,
        9.315,
        8.663,
        8.075,
        7.583,
        7.221,
        7.019,
        7.011,
        7.228,
        7.704,
        8.469,
        9.556,
        10.998,
    ]
    # fit_positions calculates also the expected position for pixel 0, not in the table
    assert fit_results.calculated_positions[1:] == pytest.approx(expected_positions, abs=1e-2)
    assert fit_results.calculated_heights[1:] == pytest.approx(expected_heights, abs=1e-2)


@pytest.fixture(scope="module")
@namedtuplefy
def data_apparent_tube_width():
    r"""Flood run to be used as input data for 'test_apparent_tube_width'"""
    return dict(
        flood_intensities=[
            [105, 96, 105, 101, 94, 102, 110, float("nan"), 105, 91],
            [110, 90, 104, 102, 99, 106, 108, float("nan"), 90, 93],
            [103, 105, 99, 101, 108, 104, 100, float("nan"), 93, 90],
            [94, 107, 102, 110, 98, 99, 101, float("nan"), 96, 109],
            [104, 101, 105, 105, 98, 110, 100, float("nan"), 109, 98],
            [101, 103, 102, 110, 106, 99, 93, float("nan"), 98, 94],
            [92, 108, float("nan"), 101, 108, 98, 105, float("nan"), 103, 98],
            [98, 92, float("nan"), 99, 101, 110, 93, float("nan"), 90, 110],
            [90, 103, 98, 104, 91, 105, 96, float("nan"), 96, 98],
            [95, 97, 109, 109, 104, 100, 95, float("nan"), 90, 97],
        ],
        wavelength_bin_boundaries=[1.0, 2.0],  # actual numbers are irrelevant
        c_tube=[
            23.6190476190476,
            23.8571428571429,
            24.5238095238095,
            24.8095238095238,
            23.9761904761905,
            24.5952380952381,
            23.8333333333333,
            float("nan"),
            23.0952380952381,
            23.2857142857143,
        ],
        c_ave=23.9550264550265,
        c_front=23.8095238095238,
        c_back=24.136902499999998,
        w_front=5.546405300938707,
        w_back=5.441993373826614,
        precision=2.0e-02,  # precision to compare reduction framework to test results
    )


@pytest.mark.parametrize(
    "workspace_with_instrument",
    [
        {
            "instrument_geometry": "n-pack",
            "n_tubes": 10,
            "n_pixels": 10,
            "diameter": 5.5e-03,
            "height": 4.2e-03,
            "spacing": 0.0,
            "x_center": 0.0,
            "y_center": 0.0,
            "z_center": 0.0,
        }
    ],
    indirect=True,
)
def test_apparent_tube_width(data_apparent_tube_width, workspace_with_instrument, temp_workspace_name):
    r"""
    Test for determining the apparent tube width, from Appendix 2, Section 2 of the master document.
    <https://www.dropbox.com/s/2mz0gy60pp9ehqm/Master%20document_110819.pdf?dl=0>

    devs - Jose Borreguero <borreguerojm@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>

    Description:
    - We use a flat detector made up of 10 tubes, each tube containing 10 pixels.
    - Even tubes make up the front panel, odd tubes make up the back panel.

    **Mantid algorithms used:**
        :ref:`DeleteWorkspaces <algm-DeleteWorkspaces-v1>`,

    **drtsans components used:**
    ~drtsans.tubecollection.TubeCollection
        <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tubecollection.py>
    """
    data = data_apparent_tube_width  # shortcut

    # Load the flood data into a Mantid workspace
    #
    flood_workspace = mtd.unique_hidden_name()  # random name for the workspace
    intensities = np.array(data.flood_intensities).reshape((10, 10, 1))
    workspace_with_instrument(
        axis_values=data.wavelength_bin_boundaries,
        intensities=intensities,
        uncertainties=np.sqrt(intensities),
        view="array",
        axis_units="wavelength",
        output_workspace=flood_workspace,
    )
    SampleLogs(flood_workspace).insert("run_number", 42)  # The flood run will have some run number, which is required
    SampleLogs(flood_workspace).insert("start_time", "2020-02-19T17:03:29.554116982")  # a start time is required
    # Calculate apparent tube widths using the pixel positions and heights of `flood_workspace` instesad of
    # retrieving them from the pixel-calibration database
    calibration = calculate_apparent_tube_width(flood_workspace, load_barscan_calibration=False)
    # Modify the pixel widths in the instrument object embedded in the workspace. We save the modifications to a
    # new workspace
    modified_flood_workspace = temp_workspace_name()  # new workspace
    calibration.apply(flood_workspace, output_workspace=modified_flood_workspace)

    # Sort the tubes assuming the ordering of the SANS instruments. In these instruments, the first half of the
    # tubes in the instrument definition file correspond to front tubes, followed by the back tubes. We
    # apply a permutation such that the order of the tubes in a SANS instrument would alternate front and
    # back tubes (hence the view='fbfb')
    collection = TubeCollection(modified_flood_workspace, "detector1").sorted(view="fbfb")
    last_tube_index = len(collection) - 1  # number of tubes, minus one because indexes begin at zero, not one
    last_pixel_index = len(collection[last_tube_index]) - 1  # number of pixels in each tube, minus one

    # compare the width of the first pixel in the first tube to the test data
    assert collection[0][0].width * 1.0e3 == pytest.approx(data.w_front, abs=data.precision)
    # compare the width of the last pixel in the last tube to the test data
    assert collection[last_tube_index][last_pixel_index].width * 1.0e3 == pytest.approx(
        data.w_back, abs=data.precision
    )

    # We do the same but now we overwrite the instrument embedded in the input workspace
    calibration = calculate_apparent_tube_width(flood_workspace, load_barscan_calibration=False)
    calibration.apply(flood_workspace)
    collection = TubeCollection(flood_workspace, "detector1").sorted(view="fbfb")
    assert collection[0][0].width * 1.0e3 == pytest.approx(data.w_front, abs=data.precision)
    last_tube_index = len(collection) - 1  # number of tubes, minus one because indexes begin at zero, not one
    last_pixel_index = len(collection[0]) - 1  # number of pixels in each tube, minus one
    assert collection[last_tube_index][last_pixel_index].width * 1.0e3 == pytest.approx(
        data.w_back, abs=data.precision
    )

    # NOTE:
    # unknown leftover workspace from workflow:
    # tubewidth_UNDEFINED_detector1_20200219:	0.0008 MB
    # manual deleting here
    DeleteWorkspace("tubewidth_UNDEFINED_detector1_20200219")


@pytest.fixture(scope="module")
@namedtuplefy
def data_generate_barscan_calibration():
    r"""Data to be used for `test_generate_barscan_calibration"""
    return dict(
        scans=[
            np.array(
                [
                    5.0,
                    96.0,
                    97.0,
                    97.0,
                    105.0,
                    21.0,
                    20.2,
                    20.0,
                    20.4,
                    104.0,
                    105.0,
                    101.0,
                    99.0,
                    98.0,
                    97.0,
                    97.0,
                    98.0,
                    103.0,
                    97.0,
                    4.0,  # tube 0
                    7.0,
                    103.0,
                    101.0,
                    104.0,
                    99.0,
                    19.2,
                    20.0,
                    20.8,
                    19.6,
                    103.0,
                    101.0,
                    102.0,
                    97.0,
                    95.0,
                    105.0,
                    105.0,
                    102.0,
                    96.0,
                    103.0,
                    7.0,  # tube 1
                    2.0,
                    97.0,
                    101.0,
                    103.0,
                    103.0,
                    19.8,
                    20.6,
                    19.8,
                    19.8,
                    98.0,
                    95.0,
                    99.0,
                    103.0,
                    100.0,
                    98.0,
                    97.0,
                    96.0,
                    95.0,
                    97.0,
                    2.0,  # tube 2
                    6.0,
                    103.0,
                    96.0,
                    98.0,
                    98.0,
                    19.6,
                    19.0,
                    20.6,
                    20.4,
                    101.0,
                    100.0,
                    97.0,
                    103.0,
                    95.0,
                    103.0,
                    101.0,
                    97.0,
                    99.0,
                    96.0,
                    8.0,
                ]
            ),  # tube3
            np.array(
                [
                    3.0,
                    97.0,
                    96.0,
                    98.0,
                    103.0,
                    99.0,
                    97.0,
                    97.0,
                    21.0,
                    19.2,
                    19.8,
                    20.6,
                    104.0,
                    99.0,
                    96.0,
                    100.0,
                    100.0,
                    95.0,
                    101.0,
                    7.0,  # tube 0
                    5.0,
                    100.0,
                    104.0,
                    102.0,
                    97.0,
                    103.0,
                    95.0,
                    104.0,
                    19.0,
                    20.2,
                    19.8,
                    19.4,
                    105.0,
                    104.0,
                    104.0,
                    95.0,
                    98.0,
                    103.0,
                    105.0,
                    4.0,  # tube 1
                    3.0,
                    95.0,
                    100.0,
                    104.0,
                    101.0,
                    95.0,
                    102.0,
                    102.0,
                    19.8,
                    20.6,
                    20.0,
                    20.6,
                    101.0,
                    103.0,
                    97.0,
                    105.0,
                    103.0,
                    96.0,
                    98.0,
                    5.0,  # tube 2
                    6.0,
                    96.0,
                    95.0,
                    97.0,
                    101.0,
                    99.0,
                    99.0,
                    102.0,
                    21.0,
                    19.2,
                    19.8,
                    19.6,
                    105.0,
                    98.0,
                    101.0,
                    97.0,
                    96.0,
                    103.0,
                    95.0,
                    5.0,
                ]
            ),  # tube3
            np.array(
                [
                    4.0,
                    97.0,
                    104.0,
                    101.0,
                    104.0,
                    101.0,
                    100.0,
                    96.0,
                    98.0,
                    95.0,
                    105.0,
                    103.0,
                    95.0,
                    19.2,
                    19.0,
                    19.4,
                    19.8,
                    104.0,
                    105.0,
                    5.0,  # tube 0
                    8.0,
                    98.0,
                    96.0,
                    104.0,
                    96.0,
                    104.0,
                    98.0,
                    96.0,
                    97.0,
                    98.0,
                    95.0,
                    99.0,
                    99.0,
                    20.8,
                    19.0,
                    19.8,
                    21.0,
                    99.0,
                    102.0,
                    4.0,  # tube 1
                    6.0,
                    104.0,
                    104.0,
                    98.0,
                    102.0,
                    105.0,
                    101.0,
                    101.0,
                    99.0,
                    101.0,
                    98.0,
                    101.0,
                    97.0,
                    19.0,
                    20.8,
                    20.4,
                    19.4,
                    95.0,
                    96.0,
                    1.0,  # tube 2
                    7.0,
                    96.0,
                    98.0,
                    104.0,
                    100.0,
                    95.0,
                    104.0,
                    100.0,
                    99.0,
                    104.0,
                    103.0,
                    103.0,
                    97.0,
                    19.4,
                    20.0,
                    19.6,
                    19.4,
                    100.0,
                    104.0,
                    4.0,
                ]
            ),
        ],  # tube3
        dcals=[
            50,
            100.0,
            150.0,
        ],  # bar positions with respect to lowest position of the bar
        # bar pixel positions allowing fitting with a fith-degree polynomial
        extended_bottom_edges=[5, 8, 10, 11, 13, 14, 16],
        # bar Y-coords allowing fitting with a fith-degree polynomial
        extended_dcals=[50, 100, 125, 135, 150, 155, 160],
        wavelength_bin_boundaries=[6.0, 7.0],
        bottom_pixels=[
            (5, 5, 5, 5),
            (8, 8, 8, 8),
            (13, 13, 13, 13),
        ],  # expected bottom pixels
        coefficients=(498.333, 27.500, -0.833, 0.000, 0.000, 0.000),
        formula="565 + {y} + 0 * {tube}",  # position of the bar in the frame of reference of the sample
        unit="mm",  # units for the pixel positions and heights
        # fitted y-coordinates for each pixel
        positions=np.array(
            [
                [
                    498.333,
                    525.0,
                    550.0,
                    573.333,
                    595.0,
                    615.0,
                    633.333,
                    650.0,
                    665.0,
                    678.333,
                    690.0,
                    700.0,
                    708.333,
                    715.0,
                    720.0,
                    723.333,
                    725.0,
                    725.0,
                    723.333,
                    720.0,
                ],
                [
                    498.333,
                    525.0,
                    550.0,
                    573.333,
                    595.0,
                    615.0,
                    633.333,
                    650.0,
                    665.0,
                    678.333,
                    690.0,
                    700.0,
                    708.333,
                    715.0,
                    720.0,
                    723.333,
                    725.0,
                    725.0,
                    723.333,
                    720.0,
                ],
                [
                    498.333,
                    525.0,
                    550.0,
                    573.333,
                    595.0,
                    615.0,
                    633.333,
                    650.0,
                    665.0,
                    678.333,
                    690.0,
                    700.0,
                    708.333,
                    715.0,
                    720.0,
                    723.333,
                    725.0,
                    725.0,
                    723.333,
                    720.0,
                ],
                [
                    498.333,
                    525.0,
                    550.0,
                    573.333,
                    595.0,
                    615.0,
                    633.333,
                    650.0,
                    665.0,
                    678.333,
                    690.0,
                    700.0,
                    708.333,
                    715.0,
                    720.0,
                    723.333,
                    725.0,
                    725.0,
                    723.333,
                    720.0,
                ],
            ]
        ),
        # fitted heights for each pixel
        heights=np.array(
            [
                [
                    27.5,
                    25.8333,
                    24.1666,
                    22.5,
                    20.8333,
                    19.1666,
                    17.5,
                    15.8333,
                    14.1666,
                    12.5,
                    10.8333,
                    9.1666,
                    7.5,
                    5.8333,
                    4.1666,
                    2.5,
                    0.8333,
                    -0.8333,
                    -2.5,
                    -4.1666,
                ],
                [
                    27.5,
                    25.8333,
                    24.1666,
                    22.5,
                    20.8333,
                    19.1666,
                    17.5,
                    15.8333,
                    14.1666,
                    12.5,
                    10.8333,
                    9.1666,
                    7.5,
                    5.8333,
                    4.1666,
                    2.5,
                    0.8333,
                    -0.8333,
                    -2.5,
                    -4.1666,
                ],
                [
                    27.5,
                    25.8333,
                    24.1666,
                    22.5,
                    20.8333,
                    19.1666,
                    17.5,
                    15.8333,
                    14.1666,
                    12.5,
                    10.8333,
                    9.1666,
                    7.5,
                    5.8333,
                    4.1666,
                    2.5,
                    0.8333,
                    -0.8333,
                    -2.5,
                    -4.1666,
                ],
                [
                    27.5,
                    25.8333,
                    24.1666,
                    22.5,
                    20.8333,
                    19.1666,
                    17.5,
                    15.8333,
                    14.1666,
                    12.5,
                    10.8333,
                    9.1666,
                    7.5,
                    5.8333,
                    4.1666,
                    2.5,
                    0.8333,
                    -0.8333,
                    -2.5,
                    -4.1666,
                ],
            ]
        ),
        precision=1.0e-03,  # one micron precision when comparing calculations with this data
    )


@pytest.mark.parametrize(
    "workspace_with_instrument",
    [
        {
            "instrument_geometry": "n-pack",
            "n_tubes": 4,
            "n_pixels": 20,
            "diameter": 5.5e-03,
            "height": 4.2e-03,
            "spacing": 0.0,
            "x_center": 0.0,
            "y_center": 0.0,
            "z_center": 0.0,
        }
    ],
    indirect=True,
)
def test_generate_barscan_calibration(
    data_generate_barscan_calibration, workspace_with_instrument, cleanfile, clean_workspace, temp_workspace_name
):
    r"""
    Test to determine pixel position and height from from Appendix 2, Section 1 of the master document.
    <https://github.com/neutrons/drtsans/blob/next/documents/Master_requirements_document.pdf>

    devs - Andrei Savici <saviciat@ornl.gov>
           Jose Borreguero <borreguerojm@ornl.gov>
    SME  - William Heller <hellerwt@ornl.gov>
           Ken Littrell <littrellkc@ornl.gov>
    """
    data = data_generate_barscan_calibration  # short nickname

    def _scan_to_file(intensities, dcal):
        r"""
        Convenience function to save the intensities of  a run holding the bar a fixed position into a Nexus file.

        Parameters
        ----------
        intensities: numpy.ndarray
            Intensity spectra for a particular run holding the bar a fixed position.
        dcal: str
            Position of the bar, in mili-meters

        Returns
        -------
        str
            Name of the Nexus file
        """
        workspace = temp_workspace_name()  # temporary workspace

        # Save intensities into a workspace endowed with an instrument. The instrument consists of a flat
        # array of four tubes, each 20 pixels.
        workspace_with_instrument(
            axis_values=data.wavelength_bin_boundaries,
            output_workspace=workspace,
            intensities=intensities.reshape(20, 4),
            view="pixel",
        )
        # Store the bar position in the workspace metadata
        AddSampleLog(
            Workspace=workspace,
            LogName="dcal_Readback",
            LogText=str(dcal),
            LogType="Number Series",
            LogUnit="mm",
        )
        # Store a fictitious run number in the workspace metadata
        SampleLogs(workspace).insert("run_number", random.randint(1, 999))
        # Store a fictitious start time for the run in the workspace metadata
        SampleLogs(workspace).insert("start_time", "2020-02-19T17:03:29.554116982")  # a start time is required
        filename = tempfile.NamedTemporaryFile("wb", suffix=".nxs").name
        cleanfile(filename)  # flag it for removal once the test finishes
        SaveNexus(InputWorkspace=workspace, Filename=filename)
        return filename

    # The plan is to have bar scans saved into nexus files
    file_names = [_scan_to_file(scan, dcal) for scan, dcal in zip(data.scans, data.dcals)]

    # Let's find the bottom-edge pixels for each tube within each run that held the bar at a fixed position,
    # and then compare with the expected data
    bottom_pixels_multi_scan = list()
    workspace = mtd.unique_hidden_name()
    # scan over each file, which contains intensities collected with the bar held at a fixed position
    for scan_index, file_name in enumerate(file_names):
        LoadNexus(Filename=file_name, OutputWorkspace=workspace)
        clean_workspace(workspace)  # delete workspace when test finishes
        # The main detector panel is called 'detector1', made up of four tubes which we gather into a TubeCollection.
        # A TubeCollection is a list of TubeSpectrum objects, each representing a physical tube. Here
        # we obtain the list of tubes for the selected double-detector-panel array.
        # The view 'decreasing X' sort the tubes by decreasing value of their corresponding X-coordinate.
        # In this view, a double detector panel looks like a single detector panel. When looking at
        # the panel standing at the sample, the leftmost tube has the highest X-coordinate, so the
        # 'decreasing X' view orders the tubes from left to right.
        collection = TubeCollection(workspace, component_name="detector1").sorted(view="decreasing X")
        # For each tube we collect the intensities (readY) and apply `find_edges` to find the bottom of the bar
        bottom_pixels = [find_edges(tube.readY.ravel()).bottom_shadow_pixel for tube in collection]
        assert bottom_pixels == pytest.approx(data.bottom_pixels[scan_index])  # compare to test data
        bottom_pixels_multi_scan.append(bottom_pixels)
    bottom_pixels_multi_scan = np.array(bottom_pixels_multi_scan)

    # Let's fit the positions of the extended bottom-edge pixels with the extended positions of the bar given in the
    # test data (data.extended_dcals), and then verify the coefficients of the fit.
    bar_formula = BarPositionFormula(formula=data.formula)
    dcals = [bar_formula.evaluate(dcal, 0) for dcal in data.extended_dcals]
    # `permissive=True` allows for unphysical fitted positions and heights, such as in this test data
    fit = fit_positions(data.extended_bottom_edges, dcals, tube_pixels=20, permissive=True)
    assert fit.coefficients == pytest.approx(data.coefficients, abs=data.precision)

    # Let's do the whole calibration. The result is a `Table` object
    calibration = calculate_barscan_calibration(
        file_names,
        component="detector1",
        order=2,
        formula=data.formula,
        permissive_fit=True,
    )
    # Compare the pixel positions and heights resulting from the calibration with the test data. A factor of 1000
    # is required because the calibration procedure assumes the input positions of the bar are mili-meters, but
    # stores the calibrated positions and heights in meters.
    assert 1000 * np.array(calibration.positions) == pytest.approx(data.positions.ravel(), abs=data.precision)
    assert 1000 * np.array(calibration.heights) == pytest.approx(data.heights.ravel(), abs=data.precision)

    # Let's use the updated pixel positions and heights to update the pixels in our temporary `workspace`
    calibration.apply(workspace)
    # Now verify the pixels positions and heights in `workspace` have indeed been updated
    collection = TubeCollection(workspace, component_name="detector1").sorted(view="decreasing X")
    positions, heights = list(), list()
    # Retrieve the pixel positions and heights for each tube of `workspace`
    for tube in collection:
        positions.append([1000 * y for y in tube.pixel_y])  # from meters to mili-meters
        heights.append([1000 * h for h in tube.pixel_heights])
    # Compare the updated pixel positions and heights of `workspace` with the test data
    assert np.array(positions) == pytest.approx(data.positions, abs=data.precision)
    assert np.array(heights) == pytest.approx(data.heights, abs=data.precision)
    # NOTE:
    # barscan_UNDEFINED_detector1_20200219:	0.00064 MB
    DeleteWorkspace("barscan_UNDEFINED_detector1_20200219")


@pytest.mark.datarepo
def test_calculate_gpsans_barscan(datarepo_dir, tmp_path, temp_workspace_name):
    r"""Calculate pixel positions and heights from a barscan, then compare to a saved barscan"""
    barscan_file_folder = path_join(datarepo_dir.gpsans, "pixel_calibration", "CG2_7465_selected_scans")
    barscan_files = [
        path_join(barscan_file_folder, filename)
        for filename in listdir(barscan_file_folder)
        if filename.startswith("nizceecbh")
    ]
    barscan_files.sort()
    calibration = calculate_barscan_calibration(barscan_files[::2])  # calibration object
    # Load save calibration and compare
    database_file = path_join(tmp_path, "saved_calibrations.json")
    calibration.save(database=database_file)
    barscan_workspace = temp_workspace_name()  # delete workspace when test finishes
    LoadNexusProcessed(barscan_files[0], OutputWorkspace=barscan_workspace)
    saved_calibration = load_calibration(
        barscan_workspace,
        "BARSCAN",
        database=database_file,
        output_workspace=temp_workspace_name(),  # table workspace
    )
    assert calibration.positions == pytest.approx(saved_calibration.positions, abs=0.1)
    assert calibration.heights == pytest.approx(saved_calibration.heights, abs=0.01)


@pytest.mark.mount_eqsans
def test_gpsans_calibration(has_sns_mount, reference_dir, temp_workspace_name):
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    # Load an events file to search a calibration for
    gpsans_file = path_join(reference_dir.gpsans, "pixel_calibration", "CG2_7465.nxs.h5")
    input_workspace = temp_workspace_name()
    LoadEventNexus(Filename=gpsans_file, OutputWorkspace=input_workspace)
    # Load a calibration
    database_file = path_join(reference_dir.gpsans, "pixel_calibration", "saved_calibrations.json")
    calibration = load_calibration(input_workspace, "BARSCAN", database=database_file)
    # Assert some data
    assert calibration.instrument == "GPSANS"
    assert calibration.positions[0:3] == pytest.approx([-0.51, -0.506, -0.502], abs=0.001)
    # Load and prepare the uncalibrated workspace
    input_workspace = mtd.unique_hidden_name()
    LoadEmptyInstrument(InstrumentName="CG2", OutputWorkspace=input_workspace)
    SampleLogs(input_workspace).insert("run_number", 8000)
    SampleLogs(input_workspace).insert("start_time", "2020-02-19T17:03:29.554116982")  # a start time is required
    calibration.apply(input_workspace)
    # Assert some data
    tube = TubeCollection(input_workspace, "detector1").sorted(view="decreasing X")[42]
    assert tube.pixel_y[0:3] == pytest.approx([-0.521, -0.517, -0.512], abs=0.001)  # units in mili-meters
    assert 1.0e3 * tube.pixel_heights[0:3] == pytest.approx([4.58, 4.57, 4.55], abs=0.01)
    # Verify as_intensities doesn't throw
    as_intensities(input_workspace)


@pytest.mark.datarepo
def test_biosans_main_detector_barscan(datarepo_dir):
    data_dir = path_join(datarepo_dir.biosans, "pixel_calibration", "runs_838_953")
    # select a subset of barscan files to perform the calibration
    first_run, last_run = 838, 953
    barscan_files = [path_join(data_dir, f"CG3_{run}.nxs") for run in range(first_run, 1 + last_run, 20)]
    calibration, addons = calculate_barscan_calibration(barscan_files, formula="{y} - 565", inspect_data=True)

    # dcal value in sample logs for the first and second workspaces
    assert addons["bar_positions"][0] == pytest.approx(1150, abs=0.001)
    assert addons["bar_positions"][1] == pytest.approx(950, abs=0.001)

    assert calibration.positions[0:3] == pytest.approx([-0.513, -0.508, -0.503], abs=0.01)


@pytest.mark.datarepo
def test_biosans_wing_detector_barscan(datarepo_dir, tmp_path, clean_workspace):
    r"""Calculate pixel positions and heights from a barscan, then compare to a saved barscan"""
    data_dir = path_join(datarepo_dir.biosans, "pixel_calibration", "runs_838_953")
    first_run, last_run = 838, 953
    detector_array = "wing_detector"  # calibration for the wing detector
    formula = "{y} - 640"  # translate from scan log value to Y-coordinate in the sample's reference frame.
    barscan_files = [path_join(data_dir, f"CG3_{run}.nxs") for run in range(first_run, 1 + last_run, 20)]
    mask_file = path_join(data_dir, "biosans_mask_bank88_tube4.xml")
    calibration = calculate_barscan_calibration(
        barscan_files, component=detector_array, formula=formula, mask=mask_file
    )
    # WARNING: this will add a small file to runtime disk, which might cause issue on the
    #          build server in the long run.
    calibration.save(database=path_join(tmp_path, "junk.json"), tablefile=path_join(tmp_path, "junk.nxs"))
    LoadNexus(barscan_files[0], OutputWorkspace="reference_workspace")
    clean_workspace("reference_workspace")
    views = calibration.as_intensities("reference_workspace")
    assert len(calibration.metadata["runnumbers"]) == len(barscan_files)
    print(views)


@pytest.mark.datarepo
def test_biosans_midrange_detector_barscan(datarepo_dir, tmp_path):
    detector_array = "midrange_detector"
    formula = "{y} - 640"  # translate from scan log value to Y-coordinate in the sample's reference frame.
    # use simulated barscan files
    barscan_file_folder = path_join(datarepo_dir.biosans, "pixel_calibration", "runs_838_953", "with_midrange")
    barscan_files = [
        path_join(barscan_file_folder, filename)
        for filename in listdir(barscan_file_folder)
        if filename.startswith("CG3_")
    ]
    barscan_files.sort()

    tube_pixels = 256
    # from visual inspection of instrument
    # mask the last tube in the midrange detector (back-midrange panel/bank104/tube4)
    closest_masked_detector_ids = list(range(106240, 106240 + tube_pixels + 1))

    # and bank 89/tube 1, the furthest from the beam
    furthest_masked_detector_ids = list(range(90112, 90112 + tube_pixels + 1))
    masked_detectors = closest_masked_detector_ids + furthest_masked_detector_ids
    # Do calibration
    calibration, addons = calculate_barscan_calibration(
        barscan_files, component=detector_array, formula=formula, mask=masked_detectors, inspect_data=True
    )
    # Assert some data
    assert calibration.instrument == "BIOSANS"
    assert calibration.positions[0:3] == pytest.approx([-0.532, -0.528, -0.524], abs=0.001)
    # Save, load, and apply calibration
    calibration.save(
        database=path_join(tmp_path, "saved_calibration.json"), tablefile=path_join(tmp_path, "saved_calibration.nxs")
    )
    calibration.as_intensities(addons["bar_workspaces"][0])
    assert len(calibration.metadata["runnumbers"]) == len(barscan_files)
    DeleteWorkspaces(addons["bar_workspaces"])


@pytest.mark.datarepo
def test_gpsans_tube_calibration(datarepo_dir, temp_workspace_name):
    r"""Calculate tube widths from a flood file"""
    flood_file = path_join(datarepo_dir.gpsans, "pixel_calibration", "CG2_8143.nxs")
    uncalibrated_workspace = temp_workspace_name()
    LoadNexus(flood_file, OutputWorkspace=uncalibrated_workspace)
    calibration = calculate_apparent_tube_width(uncalibrated_workspace, load_barscan_calibration=False)
    calibrated_workspace = temp_workspace_name()
    calibration.apply(uncalibrated_workspace, output_workspace=calibrated_workspace)

    uncalibrated_densities = _linear_density(uncalibrated_workspace)
    calibrated_densities = _linear_density(calibrated_workspace)
    assert _amplitude(calibrated_densities) / _amplitude(uncalibrated_densities) == pytest.approx(0.13, abs=0.01)
    DeleteWorkspace(calibration.table)


@pytest.mark.datarepo
def test_biosans_tube_calibration(datarepo_dir, temp_workspace_name):
    flood_file = path_join(datarepo_dir.biosans, "pixel_calibration", "flood_files", "CG3_4829.nxs")
    uncalibrated_workspace = temp_workspace_name()
    LoadNexus(flood_file, OutputWorkspace=uncalibrated_workspace)

    calibration_wing = calculate_apparent_tube_width(
        uncalibrated_workspace,
        component="wing_detector",
        load_barscan_calibration=False,
    )
    calibrated_workspace = temp_workspace_name()
    calibration_wing.apply(uncalibrated_workspace, output_workspace=calibrated_workspace)

    uncalibrated_densities = _linear_density(uncalibrated_workspace, component="wing_detector")
    calibrated_densities = _linear_density(calibrated_workspace, component="wing_detector")
    assert _amplitude(calibrated_densities) / _amplitude(uncalibrated_densities) == pytest.approx(0.38, abs=0.01)

    # cleanup
    DeleteWorkspace(calibration_wing.table)


@pytest.mark.datarepo
def test_biosans_midrange_tube_calibration(datarepo_dir, temp_workspace_name):
    flood_file = path_join(datarepo_dir.biosans, "pixel_calibration", "flood_files", "CG3_4829.nxs")
    wing_tube_calibration_workspace = LoadNexus(flood_file)

    # simulate midrange detector calibration by cloning wing detector intensities
    uncalibrated_workspace = temp_workspace_name()
    clone_component_intensities(
        wing_tube_calibration_workspace,
        output_workspace=uncalibrated_workspace,
        input_component="wing_detector",
        output_component="midrange_detector",
        copy_logs=True,
    )

    calibration_midrange = calculate_apparent_tube_width(
        uncalibrated_workspace,
        component="midrange_detector",
        load_barscan_calibration=False,
    )
    calibrated_workspace = temp_workspace_name()
    calibration_midrange.apply(uncalibrated_workspace, output_workspace=calibrated_workspace)

    uncalibrated_densities = _linear_density(uncalibrated_workspace, component="midrange_detector")
    calibrated_densities = _linear_density(calibrated_workspace, component="midrange_detector")
    assert _amplitude(calibrated_densities) / _amplitude(uncalibrated_densities) == pytest.approx(0.09, abs=0.01)

    # cleanup
    DeleteWorkspace(calibration_midrange.table)


@pytest.mark.mount_eqsans
def test_as_intensities(has_sns_mount, reference_dir):
    """Requires access to:
    /HFIR/CG2/shared/calibration/pixel_calibration.json
    """
    if not has_sns_mount:
        pytest.skip("SNS mount is not available")

    input_file = path_join(reference_dir.gpsans, "pixel_calibration", "runs_7465", "CG2_7465_55.nxs")
    LoadNexus(Filename=input_file, OutputWorkspace="scan_55")
    # "/HFIR/CG2/shared/sans-backend/data/new/ornl/sans/hfir/gpsans/pixel_calibration/runs_7465/CG2_7465_55.nxs",
    calibration = load_calibration("scan_55", "BARSCAN")
    views = calibration.as_intensities("scan_55")
    print(views)
    # NOTE:
    # mysterious leftover workspace in memory
    # barscan_GPSANS_detector1_20200103:	0.393216 MB
    # barscan_GPSANS_detector1_20200103_heights:	1.181868 MB
    # barscan_GPSANS_detector1_20200103_positions:	1.181868 MB
    # barscan_GPSANS_detector1_20200103_positions_mantid:	1.181868 MB
    # scan_55:	1.181868 MB
    DeleteWorkspace("barscan_GPSANS_detector1_20200103")
    DeleteWorkspace("barscan_GPSANS_detector1_20200103_heights")
    DeleteWorkspace("barscan_GPSANS_detector1_20200103_positions")
    DeleteWorkspace("barscan_GPSANS_detector1_20200103_positions_mantid")
    DeleteWorkspace("scan_55")


if __name__ == "__main__":
    pytest.main([__file__])
