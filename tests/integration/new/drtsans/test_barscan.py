import pytest
import numpy as np
import tempfile
import os

r""" Hyperlinks to mantid algorithms
DeleteWorkspace <https://docs.mantidproject.org/nightly/algorithms/DeleteWorkspace-v1.html>
SaveNexus <https://docs.mantidproject.org/nightly/algorithms/SaveNexus-v1.html>
AddSampleLog <https://docs.mantidproject.org/nightly/algorithms/AddSampleLog-v1.html>
"""
from mantid.simpleapi import DeleteWorkspace, SaveNexus, AddSampleLog

r"""
Hyperlinks to drtsans functions
apparent_tube_width, find_edges, fit_positions, calculate_barscan_calibration
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/barscan.py>
namedtuplefy, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
TubeCollection <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
"""  # noqa: E501
from drtsans.barscan import apparent_tube_width, find_edges, fit_positions, calculate_barscan_calibration
from drtsans.settings import namedtuplefy, unique_workspace_dundername
from drtsans.tubecollection import TubeCollection

r"""Finding the edges of the barscan in a single tube,
then calculate the position and width of the pixels

based on BarScanShadow_test_KCL_SVP.xlsx and BarScanFitTest_KCL.xlsx
Testing appendix 2.1 in the master document

devs - Andrei Savici <saviciat@ornl.gov>,
       Jose Borreguero <borreguerojm@ornl.gov>
SME  - Ken Littrell <littrellkc@ornl.gov>
"""


def test_find_edges():
    r"""Check the algorithm for finding edges in a tube as part of a barscan
    """
    # tube intensities
    intensities = np.array([2, 2, 38, 38, 38, 34, 38, 41, 35, 3, 4, 3, 3, 4,
                            30, 30, 37, 33, 31, 39, 42, 42, 2, 2, 2])
    # find edges
    edges = find_edges(intensities)
    # check expected values
    assert (edges.bottom_pixel, edges.top_pixel) == (2, 21)
    assert (edges.bottom_shadow_pixel, edges.above_shadow_pixel) == (9, 14)
    # check if the algorithm fails if there are not enough illuminated pixels
    with pytest.raises(RuntimeError, match='Faulty tube found'):
        find_edges(intensities, min_illuminated_length=len(intensities))


def test_no_shadow():
    r"""Check the failure mode for not finding shaddow
    """
    intensities = [1]*25  # all intensities are the same, no shaddow
    with pytest.raises(IndexError, match='Could not find bottom shadow edge'):
        find_edges(intensities)


def test_no_bottom_tube_edge():
    r"""Check the failure mode for not finding tube edge pixels
    """
    intensities = np.ones(25)
    intensities[14] = 3000  # all pixels but one are below the threshold
    with pytest.raises(IndexError, match='Could not find bottom tube edge'):
        find_edges(intensities)


def test_fit_positions():
    r"""Test the fitted positions and heights of pixels
    Using drtsans.barscan.fit_positions
    """
    # input edge pixels and positions
    edge_pixel = np.arange(25) + 1
    pos = [-120, -110, -98, -87, -74, -62, -49, -36, -23, -11,
           1, 12, 23, 32, 41, 49, 57, 64, 71, 78, 86, 93, 102, 110, 120]
    # fit the positions
    new_positions, new_heights = fit_positions(edge_pixel, pos, tube_pixels=26)
    # compare to the expected data
    expected_positions = [-119.659,
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
                          120.340]
    expected_heights = [9.001,
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
                        10.998]
    # fit_positions calculates also the expected position for pixel 0, not in the table
    assert new_positions[1:] == pytest.approx(expected_positions, abs=1e-2)
    assert new_heights[1:] == pytest.approx(expected_heights, abs=1e-2)


@pytest.fixture(scope='module')
@namedtuplefy
def data_apparent_tube_width():
    r"""Flood run to be used as input data for 'test_apparent_tube_width'"""
    return dict(flood_intensities=[[105, 96, 105, 101, 94, 102, 110, float('nan'), 105, 91],
                                   [110, 90, 104, 102, 99, 106, 108, float('nan'), 90, 93],
                                   [103, 105, 99, 101, 108, 104, 100, float('nan'), 93, 90],
                                   [94, 107, 102, 110, 98, 99, 101, float('nan'), 96, 109],
                                   [104, 101, 105, 105, 98, 110, 100, float('nan'), 109, 98],
                                   [101, 103, 102, 110, 106, 99, 93, float('nan'), 98, 94],
                                   [92, 108, float('nan'), 101, 108, 98, 105, float('nan'), 103, 98],
                                   [98, 92, float('nan'), 99, 101, 110, 93, float('nan'), 90, 110],
                                   [90, 103, 98, 104, 91, 105, 96, float('nan'), 96, 98],
                                   [95, 97, 109, 109, 104, 100, 95, float('nan'), 90, 97]],
                wavelength_bin_boundaries=[1.0, 2.0],  # actual numbers are irrelevant
                c_tube=[23.6190476190476, 23.8571428571429, 24.5238095238095, 24.8095238095238, 23.9761904761905,
                        24.5952380952381, 23.8333333333333, float('nan'), 23.0952380952381, 23.2857142857143],
                c_ave=23.9550264550265,
                c_front=23.8095238095238,
                c_back=24.136902499999998,
                w_front=5.46659304251795,
                w_back=5.54175869685257,
                precision=2.e-02,  # precision to compare reduction framework to test results
                )


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 10, 'n_pixels': 10,
                           'diameter': 5.5e-03, 'height': 4.2e-03, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
def test_apparent_tube_width(data_apparent_tube_width, workspace_with_instrument):
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
        <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tubecollection.py>
    """
    data = data_apparent_tube_width  # shortcut

    # Load the flood data into a Mantid workspace
    #
    flood_workspace = unique_workspace_dundername()  # random name for the workspace
    intensities = np.array(data.flood_intensities).reshape((10, 10, 1))
    workspace_with_instrument(axis_values=data.wavelength_bin_boundaries, intensities=intensities,
                              uncertainties=np.sqrt(intensities), view='array', axis_units='wavelength',
                              output_workspace=flood_workspace)

    # Find apparent tube widths and modify the pixel widths in the instrument object embedded in the workspace.
    # We save the modifications to a new workspace
    modified_flood_workspace = unique_workspace_dundername()
    apparent_tube_width(flood_workspace, output_workspace=modified_flood_workspace)

    # Sort the tubes according to the X-coordinate in decreasing value. This is the order when sitting on the
    # sample and iterating over the tubes "from left to right"
    collection = TubeCollection(modified_flood_workspace, 'detector1').sorted(view='decreasing X')
    last_tube_index = len(collection) - 1  # number of tubes, minus one because indexes begin at zero, not one
    last_pixel_index = len(collection[0]) - 1  # number of pixels in each tube, minus one

    # compare the width of the first pixel in the first tube to the test data
    assert collection[0][0].width * 1.e3 == pytest.approx(data.w_front, abs=data.precision)
    # compare the width of the last pixel in the last tube to the test data
    assert collection[last_tube_index][last_pixel_index].width * 1.e3 == pytest.approx(data.w_back, abs=data.precision)

    # We do the same but now we overwrite the instrument embedded in the input workspace
    apparent_tube_width(flood_workspace)
    collection = TubeCollection(flood_workspace, 'detector1').sorted(view='decreasing X')
    assert collection[0][0].width * 1.e3 == pytest.approx(data.w_front, abs=data.precision)
    last_tube_index = len(collection) - 1  # number of tubes, minus one because indexes begin at zero, not one
    last_pixel_index = len(collection[0]) - 1  # number of pixels in each tube, minus one
    assert collection[last_tube_index][last_pixel_index].width * 1.e3 == pytest.approx(data.w_back, abs=data.precision)

    DeleteWorkspace(modified_flood_workspace)  # flood_workspace is garbage collected upon test completion


@pytest.mark.parametrize('workspace_with_instrument',
                         [{'instrument_geometry': 'n-pack', 'n_tubes': 4, 'n_pixels': 20,
                           'diameter': 5.5e-03, 'height': 4.2e-03, 'spacing': 0.0,
                           'x_center': 0.0, 'y_center': 0.0, 'z_center': 0.0}], indirect=True)
def test_generate_calibration(workspace_with_instrument, cleanfile):
    r"""
    Test to determine pixel position and height from from Appendix 2, Section 1 of the master document.
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/documents/Master_requirements_document.pdf>

    devs - Andrei Savici <saviciat@ornl.gov>
    SME  - William Heller <hellerwt@ornl.gov>
           Ken Littrell <littrellkc@ornl.gov>
    """
    # create data files
    intensities1 = np.array([5.0, 96.0, 97.0, 97.0, 105.0, 21.0, 20.2, 20.0, 20.4, 104.0,
                             105.0, 101.0, 99.0, 98.0, 97.0, 97.0, 98.0, 103.0, 97.0, 4.0,  # tube 0
                             7.0, 103.0, 101.0, 104.0, 99.0, 19.2, 20.0, 20.8, 19.6, 103.0,
                             101.0, 102.0, 97.0, 95.0, 105.0, 105.0, 102.0, 96.0, 103.0, 7.0,  # tube 1
                             2.0, 97.0, 101.0, 103.0, 103.0, 19.8, 20.6, 19.8, 19.8, 98.0,
                             95.0, 99.0, 103.0, 100.0, 98.0, 97.0, 96.0, 95.0, 97.0, 2.0,  # tube 2
                             6.0, 103.0, 96.0, 98.0, 98.0, 19.6, 19.0, 20.6, 20.4, 101.0,
                             100.0, 97.0, 103.0, 95.0, 103.0, 101.0, 97.0, 99.0, 96.0, 8.0])  # tube3
    intensities2 = np.array([3.0, 97.0, 96.0, 98.0, 103.0, 99.0, 97.0, 97.0, 21.0, 19.2, 19.8,
                             20.6, 104.0, 99.0, 96.0, 100.0, 100.0, 95.0, 101.0, 7.0,  # tube 0
                             5.0, 100.0, 104.0, 102.0, 97.0, 103.0, 95.0, 104.0, 19.0, 20.2,
                             19.8, 19.4, 105.0, 104.0, 104.0, 95.0, 98.0, 103.0, 105.0, 4.0,  # tube 1
                             3.0, 95.0, 100.0, 104.0, 101.0, 95.0, 102.0, 102.0, 19.8, 20.6,
                             20.0, 20.6, 101.0, 103.0, 97.0, 105.0, 103.0, 96.0, 98.0, 5.0,  # tube 2
                             6.0, 96.0, 95.0, 97.0, 101.0, 99.0, 99.0, 102.0, 21.0, 19.2,
                             19.8, 19.6, 105.0, 98.0, 101.0, 97.0, 96.0, 103.0, 95.0, 5.0])  # tube3
    intensities3 = np.array([4.0, 97.0, 104.0, 101.0, 104.0, 101.0, 100.0, 96.0, 98.0, 95.0,
                             105.0, 103.0, 95.0, 19.2, 19.0, 19.4, 19.8, 104.0, 105.0, 5.0,  # tube 0
                             8.0, 98.0, 96.0, 104.0, 96.0, 104.0, 98.0, 96.0, 97.0, 98.0,
                             95.0, 99.0, 99.0, 20.8, 19.0, 19.8, 21.0, 99.0, 102.0, 4.0,  # tube 1
                             6.0, 104.0, 104.0, 98.0, 102.0, 105.0, 101.0, 101.0, 99.0,
                             101.0, 98.0, 101.0, 97.0, 19.0, 20.8, 20.4, 19.4, 95.0, 96.0, 1.0,  # tube 2
                             7.0, 96.0, 98.0, 104.0, 100.0, 95.0, 104.0, 100.0, 99.0, 104.0, 103.0,
                             103.0, 97.0, 19.4, 20.0, 19.6, 19.4, 100.0, 104.0, 4.0])  # tube3

    ws1 = workspace_with_instrument(axis_values=[6.], intensities=intensities1.reshape(20, 4), view='pixel')
    AddSampleLog(Workspace=ws1, LogName='dcal', LogText='50.', LogType='Number Series', LogUnit='mm')
    filename1 = tempfile.NamedTemporaryFile('wb', suffix='.nxs').name
    cleanfile(filename1)
    SaveNexus(InputWorkspace=ws1, Filename=filename1)
    assert (os.path.exists(filename1))

    ws2 = workspace_with_instrument(axis_values=[6.], intensities=intensities2.reshape(20, 4), view='pixel')
    AddSampleLog(Workspace=ws2, LogName='dcal', LogText='100.', LogType='Number Series', LogUnit='mm')
    filename2 = tempfile.NamedTemporaryFile('wb', suffix='.nxs').name
    cleanfile(filename2)
    SaveNexus(InputWorkspace=ws2, Filename=filename2)
    assert (os.path.exists(filename2))

    ws3 = workspace_with_instrument(axis_values=[6.], intensities=intensities3.reshape(20, 4), view='pixel')
    AddSampleLog(Workspace=ws3, LogName='dcal', LogText='150.', LogType='Number Series', LogUnit='mm')
    filename3 = tempfile.NamedTemporaryFile('wb', suffix='.nxs').name
    cleanfile(filename3)
    SaveNexus(InputWorkspace=ws3, Filename=filename3)
    assert (os.path.exists(filename3))

    barpos = calculate_barscan_calibration([filename1, filename2, filename3], 'test',
                                           order=2, formula='dcal+565', tube_pixels=20)
    assert barpos == 3
    # find edges
    edges = find_edges(intensities3[60:])
    # check expected values
    assert (edges.bottom_pixel, edges.top_pixel) == (1, 18)
    assert (edges.bottom_shadow_pixel == 5)


if __name__ == '__main__':
    pytest.main([__file__])
