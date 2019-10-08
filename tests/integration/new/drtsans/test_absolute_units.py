import pytest
from math import sqrt
import numpy as np
from mantid.simpleapi import CreateSingleValuedWorkspace, Divide
from drtsans.tof.eqsans import center_detector
from drtsans.absolute_units import empty_beam_intensity, standard_sample_scaling
from drtsans.settings import unique_workspace_name as uwn


@pytest.mark.parametrize('workspace_with_instrument',
                         [dict(name='EQSANS', Nx=17, Ny=15, dx=0.01, dy=0.01, zc=1.0)], indirect=True)
def test_empty_beam_intensity(workspace_with_instrument):
    r"""
    This test implements issue #179 and test 15B, addressing master document section 12b.

    dev - Jose Borreguero <borreguerojm@ornl.gov>
    SME - Lilin He <hel3@ornl.gov>

    **Instrument properties**:
    - 17 tubes, 15 pixels per tube.
    - pixel size is 1cm x 1cm
    - detector panel 1m away from the sample

    **Mantid algorithms used:**
    :ref:`CreateSingleValuedWorkspace <algm-CreateSingleValuedWorkspace-v1>`,
    <https://docs.mantidproject.org/nightly/algorithms/CreateSingleValuedWorkspace-v1.html>
    :ref:`Divide <algm-Divide-v1>`,
    <https://docs.mantidproject.org/nightly/algorithms/Divide-v1.html>

    **drtsans functions used:**
    ~drtsans.tof.eqsans.center_detector,
    <https://scse.ornl.gov/docs/drt/sans/drtsans/tof/eqsans/index.html>
    ~drtsans.absolute_units.empty_beam_intensity,
    <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/absolute_units.py>
    """
    # Intensities of the direct beam center on every detector pixel
    # The maximum intensity (2059.2) is located at pixel with pixel coordinates (6, 7)
    intensities = """0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.8, 0, 0, 0
0, 1.8, 5.4, 3.6, 5.4, 3.6, 5.4, 9, 14.4, 3.6, 5.4, 3.6, 3.6, 1.8, 0, 0, 1.8
3.6, 18, 21.6, 54, 91.8, 84.6, 75.6, 91.8, 72, 46.8, 45, 30.6, 14.4, 3.6, 0, 0, 0
7.2, 117, 234, 370.8, 457.2, 549, 585, 563.4, 525.6, 462.6, 441, 367.2, 259.2, 145.8, 64.8, 32.4, 21.6
84.6, 221.4, 379.8, 487.8, 473.4, 545.4, 473.4, 453.6, 345.6, 331.2, 280.8, 205.2, 131.4, 59.4, 19.8, 5.4, 1.8
174.6, 657, 930.6, 1373.4, 1459.8, 1411.2, 1405.8, 1339.2, 1159.2, 1035, 950.4, 748.8, 487.8, 342, 169.2, 41.4, 18
82.8, 333, 487.8, 576, 925.2, 837, 907.2, 739.8, 817.2, 720, 563.4, 545.4, 379.8, 275.4, 147.6, 61.2, 23.4
127.8, 583.2, 1040.4, 1396.8, 1702.8, 1956.6, 2059.2, 1929.6, 1841.4, 1825.2, 1582.2, 1281.6, 1078.2, 658.8, 378, 185.4, 79.2
228.6, 567, 757.8, 923.4, 1065.6, 1006.2, 927, 761.4, 738, 590.4, 529.2, 387, 217.8, 136.8, 66.6, 10.8, 7.2
257.4, 734.4, 1096.2, 1450.8, 1542.6, 1629, 1587.6, 1427.4, 1305, 1054.8, 939.6, 730.8, 516.6, 288, 111.6, 48.6, 21.6
63, 221.4, 358.2, 495, 549, 608.4, 626.4, 547.2, 473.4, 424.8, 311.4, 284.4, 162, 81, 36, 28.8, 10.8
43.2, 214.2, 325.8, 478.8, 549, 552.6, 570.6, 586.8, 457.2, 403.2, 363.6, 207, 201.6, 82.8, 27, 18, 5.4
10.8, 19.8, 41.4, 70.2, 115.2, 95.4, 106.2, 79.2, 88.2, 55.8, 41.4, 18, 9, 3.6, 1.8, 1.8, 0
0, 1.8, 9, 14.4, 16.2, 25.2, 19.8, 25.2, 25.2, 10.8, 3.6, 5.4, 1.8, 0, 0, 0, 0
0, 0, 0, 0, 0, 0, 0, 1.8, 0, 0, 0, 0, 0, 0, 0, 0, 0"""  # noqa: E501
    center_xy = (6, 7)  # center of the bin, in pixel coordinates
    pixels_per_tube = 15  # number of pixels in any given tube

    # save the intensities in a numpy.ndarray.
    intensities = [float(x) for line in intensities.split('\n') for x in line.split(',')]
    intensities = np.array(intensities, dtype=float).reshape((15, 17, 1))

    # Create a Mantid workspace with an embedded instrument and store the attenuated empty beam intensities
    wavelength_bin = [2.5, 3.0]  # some arbitrary wavelength bin
    empty_beam_workspace = workspace_with_instrument(axis_values=wavelength_bin, intensities=intensities, view='array')

    # Confirm the intensity of the center beam (pixel with pixel coordinates (6, 7))
    center_pixel_id = center_xy[0] * pixels_per_tube + center_xy[1]
    maximum_intensity = empty_beam_workspace.readY(center_pixel_id)[0]
    assert maximum_intensity == pytest.approx(2059.2)

    # Find the coordinates of the center pixel, then translate the detector such that the coordinates of the
    # center pixel lie now at the origin of the XY plane
    spectrum_info = empty_beam_workspace.spectrumInfo()
    x, y, z = spectrum_info.position(center_pixel_id)
    center_detector(empty_beam_workspace, center_x=x, center_y=y)

    # Check that the center pixel is now at the origin in the XY plane
    assert spectrum_info.position(center_pixel_id) == pytest.approx((0.0, 0.0, 1.0))

    # Find the intensity of the attenuated empty beam within a radius of 60mm around the beam center,
    # and remove the attenuation with the input attenuation coefficient
    empty_beam_intensity(empty_beam_workspace, beam_radius=60, unit='mm', output_workspace='empty_beam_intensity',
                         attenuator_coefficient=1./30, attenuator_error=0.01/30)

    # Store the input sample intensity in a Mantid workspace
    CreateSingleValuedWorkspace(DataValue=100000, ErrorValue=sqrt(100000), OutputWorkspace='sample_intensity')

    # The normalized sample intensity is obtain by dividing with the empty beam intensity
    normalized_intensity = Divide(LHSWorkspace='sample_intensity', RHSWorkspace='empty_beam_intensity')

    # Compare with the values written in the test
    assert normalized_intensity.readY(0)[0] == pytest.approx(0.04229, abs=1e-05)
    assert normalized_intensity.readE(0)[0] == pytest.approx(0.00046, abs=1e-05)


def test_standard_sample_measurement():
    '''
        Tests normalization with a calibrated standard sample as described in the master document
        section 12.3
        dev - Steven Hahn <hahnse@ornl.gov>
        SME - Changwoo Do <doc1@ornl.gov>
    '''
    F_std = 450.  # value from supplied test
    F_std_err = 10.  # value from supplied test
    F_std_ws = CreateSingleValuedWorkspace(DataValue=F_std, ErrorValue=F_std_err, OutputWorkspace=uwn())
    F = 10.  # value from supplied test
    F_err = 2.  # value from supplied test
    F_ws = CreateSingleValuedWorkspace(DataValue=F, ErrorValue=F_err, OutputWorkspace=uwn())
    Iq = 100.  # value from supplied test
    Iq_err = np.sqrt(Iq)
    Iq_ws = CreateSingleValuedWorkspace(DataValue=Iq, ErrorValue=Iq_err, OutputWorkspace=uwn())
    # perform calculation done by function standard_sample_scaling
    Iq_abs = Iq / F * F_std
    # calculate uncertainty as described in the supplied test. Symbolically this is identical to the calculation below,
    # but, due to numerical error, differs both the calculation below and Mantid's result by ~10%
    # err1 = F_std**2 / F**4 * F_err**2 * Iq**2
    # err2 = F_std_err**2 / F**2 * Iq**2
    # err3 = F_std**2 / F**2 * Iq_err*2
    # Iq_abs_err = np.sqrt(err1 + err2 + err3)
    err1 = F_err**2 / F**2
    err2 = F_std_err**2 / F_std**2
    err3 = Iq_err**2 / Iq**2
    Iq_abs_err = Iq / F * F_std * np.sqrt(err1 + err2 + err3)
    # run absolute_units.standard_sample_scaling
    Iq_abs_ws = standard_sample_scaling(Iq_ws, F_ws, F_std_ws)
    # check results
    assert Iq_ws.dataY(0)[0] == pytest.approx(Iq)
    assert Iq_abs_ws.dataY(0)[0] == pytest.approx(Iq_abs)
    assert Iq_abs_ws.dataE(0)[0] == pytest.approx(Iq_abs_err)
