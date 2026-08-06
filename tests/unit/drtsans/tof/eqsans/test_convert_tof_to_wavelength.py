from dataclasses import replace
import os

import pytest
from mantid.simpleapi import LoadEventNexus

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py
from drtsans.samplelogs import SampleLogs

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/correct_frame.py
from drtsans.tof.eqsans.correct_frame import (
    convert_to_wavelength,
    transform_to_wavelength,
    IncompatibleWavelengthBandsError,
    MONOCHROMATIC_PV,
    TransmittedBands,
    WAVELENGTH_BAND_DIFF_TOLERANCE,
)
from drtsans.wavelength import Wband, from_tof


def add_frame_skipping_log(ws):
    samplelog = SampleLogs(ws)
    samplelog.insert("is_frame_skipping", False)


def add_monochromatic_log(ws, value=1):
    """Emulate the data acquisition system flagging the run as taken in monochromatic mode"""
    SampleLogs(ws).insert(MONOCHROMATIC_PV, value)


@pytest.mark.parametrize(
    "generic_workspace",
    [
        {
            "dx": 0.005,
            "dy": 0.004,
            "zc": 5.0,
            "l1": 14.0,
            "axis_units": "tof",
            "axis_values": [15432.0],
        }
    ],
    indirect=True,
)
def test_william(generic_workspace, clean_workspace):
    """Test the conversion of time-of-flight to wavelength
    in master document section 3.3
    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - William Heller <hellerwt@ornl.gov>"""
    # set up input workspace in the first frame
    ws = generic_workspace  # friendly name
    clean_workspace(ws)
    add_frame_skipping_log(ws)

    # run the drt-sans version
    ws = convert_to_wavelength(input_workspace=generic_workspace)

    # make sure the unit is wavelength
    assert ws.getAxis(0).getUnit().caption() == "Wavelength"

    # get information for detector pixel positions
    specInfo = ws.spectrumInfo()
    source_sample = specInfo.l1()  # in meters

    # verify the individual wavelength values
    for i in range(4):
        # distance to detector pixel in meters
        sample_detector = specInfo.l2(i)
        # Calculate expected wavelength using drtsans.wavelength.from_tof
        expected_wavelength = from_tof(15432.0, distance=source_sample + sample_detector)

        # With FullBinsOnly=True, bin edges are adjusted for complete bins
        # Allow small tolerance for bin grid adjustments
        assert ws.dataX(i)[0] == pytest.approx(expected_wavelength, rel=0.02)


TOF = [12345.0, 12346.0]


@pytest.mark.parametrize(
    "generic_workspace",
    [
        {
            "dx": 0.005,
            "dy": 0.004,
            "zc": 2.5,
            "l1": 10.1,
            "axis_units": "tof",
            "axis_values": TOF,
        }
    ],
    indirect=True,
)
def test_shuo(generic_workspace, clean_workspace):
    """Test the conversion of time-of-flight to wavelength
    in master document section 3.3
    dev - Pete Peterson <petersonpf@ornl.gov>
    SME - Shuo Qian"""
    # set up input workspace in the first frame
    ws = generic_workspace  # friendly name
    clean_workspace(ws)
    add_frame_skipping_log(ws)

    # run the drt-sans version
    ws = convert_to_wavelength(input_workspace=generic_workspace)

    # make sure the unit is wavelength
    assert ws.getAxis(0).getUnit().caption() == "Wavelength"

    # get information for detector pixel positions
    specInfo = ws.spectrumInfo()
    source_sample = specInfo.l1()  # in meters

    # verify the individual wavelength values
    for i in range(4):
        # distance to detector pixel in meters
        sample_detector = specInfo.l2(i)
        # Calculate expected wavelength for first bin edge using drtsans.wavelength.from_tof
        expected_wavelength_0 = from_tof(TOF[0], distance=source_sample + sample_detector)

        # With FullBinsOnly=True, bin edges are adjusted for complete bins
        # Check the first bin edge matches the expected conversion from TOF
        assert ws.dataX(i)[0] == pytest.approx(expected_wavelength_0, rel=0.02)

        # Verify wavelength values are in reasonable range (monotonically increasing)
        assert ws.dataX(i)[0] < ws.dataX(i)[1]


# Workspace spanning, after conversion to wavelength, a range wider than `MONOCHROMATIC_BAND`
MONOCHROMATIC_WORKSPACE = {
    "dx": 0.005,
    "dy": 0.004,
    "zc": 5.0,
    "l1": 14.0,
    "axis_units": "tof",
    "axis_values": [10000.0, 20000.0],
}

# Wavelength band transmitted by the choppers, in Angstroms
MONOCHROMATIC_BAND = Wband(2.5, 3.5)


@pytest.mark.parametrize("generic_workspace", [MONOCHROMATIC_WORKSPACE], indirect=True)
def test_convert_to_wavelength_monochromatic(generic_workspace, clean_workspace):
    """In monochromatic mode the requested bin width is overridden by a single bin spanning
    the whole transmitted band"""
    ws = generic_workspace
    clean_workspace(ws)
    add_frame_skipping_log(ws)
    add_monochromatic_log(ws)

    bands = TransmittedBands(lead=MONOCHROMATIC_BAND, skip=None)
    ws = convert_to_wavelength(ws, bands=bands, bin_width=0.1)

    assert ws.dataX(0) == pytest.approx([MONOCHROMATIC_BAND.min, MONOCHROMATIC_BAND.max])
    band_width = MONOCHROMATIC_BAND.max - MONOCHROMATIC_BAND.min  # Angstrom
    assert SampleLogs(ws).single_value("wavelength_bin_width") == pytest.approx(band_width)


@pytest.mark.parametrize("generic_workspace", [MONOCHROMATIC_WORKSPACE], indirect=True)
def test_convert_to_wavelength_not_monochromatic(generic_workspace, clean_workspace):
    """Without the MCON16 sample log, the requested bin width is honored"""
    ws = generic_workspace
    clean_workspace(ws)
    add_frame_skipping_log(ws)

    bands = TransmittedBands(lead=MONOCHROMATIC_BAND, skip=None)
    ws = convert_to_wavelength(ws, bands=bands, bin_width=0.1)

    wavelength_bin_edges = ws.dataX(0)
    assert len(wavelength_bin_edges) > 2  # more than a single bin
    assert wavelength_bin_edges[1] - wavelength_bin_edges[0] == pytest.approx(0.1)


@pytest.mark.parametrize("generic_workspace", [MONOCHROMATIC_WORKSPACE], indirect=True)
def test_convert_to_wavelength_monochromatic_frame_skipping(generic_workspace, clean_workspace):
    """Monochromatic mode and frame-skipping mode are mutually exclusive chopper settings"""
    ws = generic_workspace
    clean_workspace(ws)
    SampleLogs(ws).insert("is_frame_skipping", True)
    add_monochromatic_log(ws)

    bands = TransmittedBands(lead=Wband(2.5, 3.0), skip=Wband(3.2, 3.5))
    with pytest.raises(ValueError, match="incompatible with frame-skipping"):
        convert_to_wavelength(ws, bands=bands)


@pytest.mark.parametrize("generic_workspace", [MONOCHROMATIC_WORKSPACE], indirect=True)
def test_convert_to_wavelength_monochromatic_unknown_bands(generic_workspace, clean_workspace):
    """The single bin cannot be constructed if the transmitted band is neither passed nor logged"""
    ws = generic_workspace
    clean_workspace(ws)
    add_frame_skipping_log(ws)
    add_monochromatic_log(ws)

    with pytest.raises(ValueError, match="requires the wavelength bands"):
        convert_to_wavelength(ws)


@pytest.mark.datarepo
def test_transform_to_wavelength_compatible_bands(datarepo_dir, temp_workspace_name):
    file_name = os.path.join(datarepo_dir.eqsans, "EQSANS_86217.nxs.h5")
    ws_tof = LoadEventNexus(file_name, OutputWorkspace=temp_workspace_name())
    _, bands = transform_to_wavelength(ws_tof)

    small_band_diff = 0.1 * WAVELENGTH_BAND_DIFF_TOLERANCE

    # create band with a lead max that is different within the tolerance
    new_lead = Wband(bands.lead.min, bands.lead.max + small_band_diff)
    new_bands_lead = replace(bands, lead=new_lead)

    _, bands_out = transform_to_wavelength(ws_tof, bands=new_bands_lead)
    assert bands_out.almost_equal(new_bands_lead)

    # create band with a skip min that is different within the tolerance
    new_skip = Wband(bands.skip.min + small_band_diff, bands.skip.max)
    new_bands_skip = replace(bands, skip=new_skip)

    _, bands_out = transform_to_wavelength(ws_tof, bands=new_bands_skip)
    assert bands_out.almost_equal(new_bands_skip)


@pytest.mark.datarepo
def test_transform_to_wavelength_raises_error_for_incompatible_bands(datarepo_dir, temp_workspace_name):
    file_name = os.path.join(datarepo_dir.eqsans, "EQSANS_86217.nxs.h5")
    ws_tof = LoadEventNexus(file_name, OutputWorkspace=temp_workspace_name())
    _, bands = transform_to_wavelength(ws_tof)

    # create band with a lead max that is incompatible with the lead max in the data
    incompatible_lead = Wband(2.480, 6.50)
    incompatible_bands_lead = replace(bands, lead=incompatible_lead)

    with pytest.raises(IncompatibleWavelengthBandsError):
        transform_to_wavelength(ws_tof, bands=incompatible_bands_lead)

    # create band with a skip min that is incompatible with the skip min in the data
    incompatible_skip = Wband(10.4, 15.230)
    incompatible_bands_skip = replace(bands, skip=incompatible_skip)

    with pytest.raises(IncompatibleWavelengthBandsError):
        transform_to_wavelength(ws_tof, bands=incompatible_bands_skip)


if __name__ == "__main__":
    pytest.main([__file__])
