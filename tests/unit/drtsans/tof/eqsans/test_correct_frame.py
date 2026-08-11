from collections import namedtuple
from os.path import join as pjoin

from mantid.kernel import DateAndTime, amend_config, FloatTimeSeriesProperty
from mantid.simpleapi import AddSampleLog, Load, CreateWorkspace, CreateSampleWorkspace
import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose
import pytest
from pytest import approx

from drtsans import wavelength as sans_wavelength
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import correct_frame
from drtsans.geometry import source_detector_distance

BandsTuple = namedtuple("BandsTuple", "lead skip")


@pytest.mark.datarepo
@pytest.mark.parametrize(
    "filename, lead_range, skip_range",
    [
        # lead_range and skip_range are (minimum, maximum) wavelengths in Angstrom
        # four chopper configuration (before 2026)
        ("EQSANS_101595.nxs.h5", (1.95, 6.08), None),
        ("EQSANS_86217.nxs.h5", (2.45, 6.71), (10.96, 15.16)),  # frame skipping mode
        # six chopper configuration, offsets effective 2026-03-04 onward
        # (no test data currently covers the 2026-01-01..2026-03-03 sub-era)
        ("EQSANS_176973.nxs.h5", (11.95, 14.89), None),
        ("EQSANS_176937.nxs.h5", (2.45, 6.04), None),
        ("EQSANS_178264.nxs.h5", (2.45, 6.06), (9.66, 13.32)),  # frame skipping mode
    ],
)
def test_transmitted_bands(datarepo_dir, clean_workspace, filename, lead_range, skip_range):
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename=filename)
        clean_workspace(ws)
        bands = correct_frame.transmitted_bands(ws)
        assert_almost_equal((bands.lead.min, bands.lead.max), lead_range, decimal=2)
        if skip_range is not None:
            assert_almost_equal((bands.skip.min, bands.skip.max), skip_range, decimal=2)
        else:
            assert bands.skip is None


@pytest.mark.datarepo
def test_transmitted_bands_zero_speed_choppers(datarepo_dir, clean_workspace):
    """Test that the transmitted bands are correctly calculated when some choppers have zero speed."""
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename="EQSANS_86217.nxs.h5")
        AddSampleLog(ws, "Speed5", "0", LogType="Number Series")
        AddSampleLog(ws, "Phase5", "0", LogType="Number Series")
        AddSampleLog(ws, "Speed6", "0", LogType="Number Series")
        AddSampleLog(ws, "Phase6", "0", LogType="Number Series")
        # overwrite start_time log to simulate run with new chopper configuration
        AddSampleLog(ws, "start_time", "2026-03-05T05:49:47.754251666", LogType="String")
        clean_workspace(ws)

        run = ws.mutableRun()
        # To simulate the new chopper configuration, we need to update the "Phase" logs of the first 4 choppers
        # to match the phase offsets in the new configuration
        # (the offset differences are the difference between the old and new offsets for frame-skip mode)
        offset_difference = [10884.92, 10790.8, 9738.12, 9771.38]
        for i in range(4):
            phase_log_name = "Phase{}".format(i + 1)
            phase_log = run.getProperty(phase_log_name)
            times = phase_log.times
            values = phase_log.value
            new_phase_log = FloatTimeSeriesProperty(phase_log_name)
            for t, v in zip(times, values):
                new_phase_log.addValue(t, v + offset_difference[i])
            run.addProperty(phase_log_name, new_phase_log, True)

        bands = correct_frame.transmitted_bands(ws)
        assert_almost_equal((bands.lead.min, bands.lead.max), (2.45, 6.73), decimal=2)
        assert_almost_equal((bands.skip.min, bands.skip.max), (11.01, 15.21), decimal=2)


@pytest.mark.datarepo
def test_transmitted_bands_clipped(datarepo_dir, clean_workspace):
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename="EQSANS_86217.nxs.h5")
        clean_workspace(ws)
        sdd = source_detector_distance(ws, unit="m")
        bands_0 = correct_frame.transmitted_bands_clipped(ws, sdd, 0.0, 0.0)
        lwc, hwc = (0.139, 0.560)  # expected clippings
        # With no interior clipping
        bands = correct_frame.transmitted_bands_clipped(ws, sdd, 500, 2000, interior_clip=False)
        # Check clippings for the lead pulse
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0), 0.01)
        # Check clippings for the skip pulse
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0, b2_0 - hwc), 0.01)
        # With interior clipping
        bands = correct_frame.transmitted_bands_clipped(ws, sdd, 500, 2000, interior_clip=True)
        b1_0, b2_0 = bands_0.lead.min, bands_0.lead.max
        b1, b2 = bands.lead.min, bands.lead.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)
        b1_0, b2_0 = bands_0.skip.min, bands_0.skip.max
        b1, b2 = bands.skip.min, bands.skip.max
        assert (b1, b2) == approx((b1_0 + lwc, b2_0 - hwc), 0.01)


def test_is_monochromatic(temp_workspace_name, clean_workspace):
    """It is the value of the MCON16 log, not its presence, that flags monochromatic mode"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    # runs predating the process variable carry no such log
    assert correct_frame.is_monochromatic(ws) is False
    SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_PV, 0)
    assert correct_frame.is_monochromatic(ws) is False
    SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_PV, 1)
    assert correct_frame.is_monochromatic(ws) is True


def insert_monochromatic_band_logs(ws, center, spread):
    """Emulate the DAS recording the wavelength band requested of the chopper controller

    Parameters
    ----------
    center: float
        middle of the band, in Angstrom
    spread: float
        width of the band, as a percent of `center`
    """
    SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_CENTER_PV, center)
    SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_SPREAD_PV, spread)


@pytest.mark.parametrize(
    "center, spread, expected",
    [
        (2.5, 15.0, (2.3125, 2.6875)),
        (1.0, 3.0, (0.985, 1.015)),
    ],
)
def test_band_from_logs(temp_workspace_name, clean_workspace, center, spread, expected):
    """The requested band spans center * (1 -+ spread / 200)"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    insert_monochromatic_band_logs(ws, center, spread)
    band = correct_frame.band_from_logs(ws)
    assert (band.min, band.max) == approx(expected, abs=1.0e-9)


@pytest.mark.parametrize("present", [None, "center", "spread"])
def test_band_from_logs_missing_logs(temp_workspace_name, clean_workspace, present):
    """Runs predating the two process variables carry neither, even when flagged monochromatic"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    if present == "center":
        SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_CENTER_PV, 2.5)
    elif present == "spread":
        SampleLogs(ws).insert(correct_frame.MONOCHROMATIC_SPREAD_PV, 15.0)

    with pytest.raises(correct_frame.MissingMonochromaticLogs, match="not found"):
        correct_frame.band_from_logs(ws)


@pytest.mark.parametrize(
    "center, spread",
    [
        (2.5, 0.0),  # zero width, would divide by zero when computing the overlap
        (2.5, -3.0),  # lower boundary above the upper one
        (2.5, 250.0),  # negative lower boundary
        (float("nan"), 10.0),
        (2.5, float("nan")),
        (2.5, float("inf")),
    ],
)
def test_band_from_logs_degenerate(temp_workspace_name, clean_workspace, center, spread):
    """Log values that do not describe a usable band are rejected with a clear error"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    insert_monochromatic_band_logs(ws, center, spread)

    with pytest.raises(correct_frame.DegenerateMonochromaticBand, match="unusable band"):
        correct_frame.band_from_logs(ws)


def test_band_from_logs_full_spread(temp_workspace_name, clean_workspace):
    """A spread of exactly 200 percent spans zero to twice the center, which is still usable"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    insert_monochromatic_band_logs(ws, 2.5, 200.0)

    band = correct_frame.band_from_logs(ws)

    assert (band.min, band.max) == approx((0.0, 5.0), abs=1.0e-9)


# Requested band for a 2.5 Angstrom center and a 10% spread, hence 2.375-2.625 Angstrom and
# 0.25 Angstrom wide. The geometric bands below cover a known fraction of it.
@pytest.mark.parametrize(
    "phased_for, expected_overlap, severity",
    [
        (sans_wavelength.Wband(2.375, 2.625), 1.00, "information"),  # exact match
        (sans_wavelength.Wband(2.3, 2.7), 1.00, "information"),  # requested band fully inside
        (sans_wavelength.Wband(2.4125, 2.625), 0.85, "warning"),  # 0.2125 of 0.25
        (sans_wavelength.Wband(2.5, 2.625), 0.50, "error"),  # 0.125 of 0.25
        (sans_wavelength.Wband(3.421, 3.562), 0.00, "error"),  # disjoint, as in run 186249
    ],
)
def test_verify_monochromatic_band(
    temp_workspace_name, clean_workspace, mocker, phased_for, expected_overlap, severity
):
    """The overlap is the fraction of the requested band the choppers are phased for, and it is
    reported at a severity increasing with the disagreement"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    insert_monochromatic_band_logs(ws, 2.5, 10.0)
    mock_logger = mocker.patch("drtsans.tof.eqsans.correct_frame.logger")

    if severity == "error":
        with pytest.raises(ValueError, match="of the requested"):
            correct_frame.verify_monochromatic_band(ws, phased_for)
    else:
        overlap = correct_frame.verify_monochromatic_band(ws, phased_for)
        assert overlap == approx(expected_overlap, abs=1.0e-3)

    # the severity itself is part of the contract, not just whether the call raised
    assert getattr(mock_logger, severity).call_count == 1
    for other in {"information", "warning", "error"} - {severity}:
        assert getattr(mock_logger, other).call_count == 0


def test_verify_monochromatic_band_unverifiable(temp_workspace_name, clean_workspace, mocker):
    """Without the two process variables the check is skipped with a warning, not an exception"""
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    mock_logger = mocker.patch("drtsans.tof.eqsans.correct_frame.logger")

    assert correct_frame.verify_monochromatic_band(ws, sans_wavelength.Wband(2.375, 2.625)) is None

    assert mock_logger.warning.call_count == 1
    assert "cannot verify" in mock_logger.warning.call_args.args[0]
    assert mock_logger.error.call_count == 0


def test_verify_monochromatic_band_degenerate(temp_workspace_name, clean_workspace, mocker):
    """A zero-width requested band skips the check with a warning rather than failing the run.

    Computing the overlap would divide by the zero width of the requested band.
    """
    ws = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    clean_workspace(ws)
    insert_monochromatic_band_logs(ws, 2.5, 0.0)
    mock_logger = mocker.patch("drtsans.tof.eqsans.correct_frame.logger")

    assert correct_frame.verify_monochromatic_band(ws, sans_wavelength.Wband(2.375, 2.625)) is None

    assert mock_logger.warning.call_count == 1
    assert "cannot verify" in mock_logger.warning.call_args.args[0]
    # not an error, and above all not an exception
    assert mock_logger.error.call_count == 0


@pytest.mark.datarepo
def test_transform_to_wavelength_clips_polychromatic(datarepo_dir, clean_workspace):
    """In the usual polychromatic mode the TOF clippings trim both edges of the band"""
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename="EQSANS_176937.nxs.h5")  # no MCON16 sample log
        clean_workspace(ws)
        assert correct_frame.is_monochromatic(ws) is False

        unclipped = correct_frame.transmitted_bands(ws)
        sdd = source_detector_distance(ws, unit="m")
        low_clip = sans_wavelength.from_tof(500.0, distance=sdd)  # Angstrom
        high_clip = sans_wavelength.from_tof(2000.0, distance=sdd)  # Angstrom

        _, bands = correct_frame.transform_to_wavelength(ws, low_tof_clip=500, high_tof_clip=2000)

        assert bands.lead.min == approx(unclipped.lead.min + low_clip, abs=1.0e-4)
        assert bands.lead.max == approx(unclipped.lead.max - high_clip, abs=1.0e-4)
        # the clippings actually used are recorded in the logs
        sample_logs = SampleLogs(ws)
        assert sample_logs.single_value("low_tof_clip") == approx(500.0)
        assert sample_logs.single_value("high_tof_clip") == approx(2000.0)


@pytest.mark.datarepo
def test_transform_to_wavelength_monochromatic_ignores_clips(datarepo_dir, clean_workspace):
    """In monochromatic mode the TOF clippings are discarded, preserving the narrow band.

    Clipping this run by the schema defaults of 500 and 2000 micro seconds would leave only
    9.59-9.85 Angstrom, just over a quarter of the transmitted band.
    """
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename="EQSANS_177103.nxs.h5")
        clean_workspace(ws)
        assert correct_frame.is_monochromatic(ws) is True

        _, bands = correct_frame.transform_to_wavelength(ws, low_tof_clip=500, high_tof_clip=2000)

        assert_almost_equal((bands.lead.min, bands.lead.max), (9.446, 10.410), decimal=2)
        assert bands.skip is None
        # the clipped band would be 0.26 Angstrom wide
        assert bands.lead.max - bands.lead.min > 0.9
        # the choppers are phased for a round 9.5-10.5 Angstrom; the emission delay shifts the
        # band that reaches the sample towards shorter wavelengths
        assert_almost_equal(
            (correct_frame.geometric_band(ws).min, correct_frame.geometric_band(ws).max),
            (9.500, 10.500),
            decimal=2,
        )
        # the override is recorded in the logs, whatever the configuration asked for
        sample_logs = SampleLogs(ws)
        assert sample_logs.single_value("low_tof_clip") == approx(0.0)
        assert sample_logs.single_value("high_tof_clip") == approx(0.0)


@pytest.mark.datarepo
def test_transform_to_wavelength_monochromatic_unverifiable(datarepo_dir, clean_workspace, mocker):
    """A monochromatic run predating MCWL16 reduces normally, with the band check skipped.

    This run is flagged monochromatic but carries neither MCWL16 nor MCWLSpread16, so there is
    no requested band to compare the chopper settings against.
    """
    with amend_config(data_dir=datarepo_dir.eqsans):
        ws = Load(Filename="EQSANS_177103.nxs.h5")
        clean_workspace(ws)
        assert correct_frame.is_monochromatic(ws) is True
        with pytest.raises(correct_frame.MissingMonochromaticLogs):
            correct_frame.band_from_logs(ws)

        mock_logger = mocker.patch("drtsans.tof.eqsans.correct_frame.logger")
        _, bands = correct_frame.transform_to_wavelength(ws, low_tof_clip=500, high_tof_clip=2000)

        assert_almost_equal((bands.lead.min, bands.lead.max), (9.446, 10.410), decimal=2)
        warnings = [call.args[0] for call in mock_logger.warning.call_args_list]
        assert any("cannot verify" in message for message in warnings), f"no skip warning in {warnings}"


@pytest.mark.datarepo
def test_log_tof_structure(datarepo_dir, temp_workspace_name):
    # reuse the same file
    file_name = pjoin(datarepo_dir.eqsans, "test_chopper", "EQSANS_92353_no_events.nxs")
    for ny, refv in ((False, 30833), (True, 28333)):
        ws = Load(file_name, OutputWorkspace=temp_workspace_name())
        correct_frame.log_tof_structure(ws, 500, 2000, interior_clip=ny)
        sl = SampleLogs(ws)
        assert sl.tof_frame_width.value == approx(33333, abs=1.0)
        assert sl.tof_frame_width_clipped.value == approx(refv, abs=1.0)


def test_band_structure_logs(temp_workspace_name):
    w = CreateWorkspace([0], [0], OutputWorkspace=temp_workspace_name())
    with pytest.raises(RuntimeError, match="Band structure not found in the logs"):
        correct_frame.metadata_bands(w)
    SampleLogs(w).insert("is_frame_skipping", 0)
    correct_frame.log_band_structure(w, BandsTuple(sans_wavelength.Wband(1.5, 2.42), None))
    bands = correct_frame.metadata_bands(w)
    assert bands.lead.min, bands.lead.max == approx(1.5, 2.42)
    SampleLogs(w).insert("is_frame_skipping", 1)
    with pytest.raises(RuntimeError, match="Bands from the skipped pulse missing in the logs"):
        correct_frame.metadata_bands(w)
    correct_frame.log_band_structure(
        w,
        BandsTuple(sans_wavelength.Wband(1.5, 2.42), sans_wavelength.Wband(6.07, 10.01)),
    )
    bands = correct_frame.metadata_bands(w)
    assert bands.lead.min, bands.lead.max == approx(1.5, 2.42)
    assert bands.skip.min, bands.skip.max == approx(6.07, 10.01)


def test_correct_emission_time_60Hz(clean_workspace):
    # excepted wavelengths
    expected_wl = np.arange(1.05, 5.30, 0.1)
    starting_tof = [
        4184.43261874537,
        4585.6620961358,
        4985.12109583935,
        5383.33670855604,
        5780.65753233583,
        6177.18810337877,
        6572.86455843483,
        6967.671528804,
        7362.0002659363,
        7757.14799863174,
        8146.18938057224,
        8537.85660707455,
        8929.51955517686,
        9321.17826927916,
        9712.83279378146,
        10104.4831730838,
        10496.1294515861,
        10887.7716736884,
        11279.4098837907,
        11671.044126293,
        12062.6744455953,
        12454.3008860976,
        12845.9234921999,
        13237.5423083022,
        13629.1573788045,
        14020.7687481068,
        14412.3764606091,
        14803.9805607114,
        15195.5810928137,
        15587.178101316,
        15978.7716306183,
        16370.3617251206,
        16761.9484292229,
        17153.5317873253,
        17545.1118438275,
        17936.6886431299,
        18328.2622296322,
        18719.8326477345,
        19111.3999418368,
        19502.9641563391,
        19894.5253356414,
        20286.0835241437,
        20677.638766246,
    ]

    # Make a simple workspace with correct distances and add tofs to it
    w = CreateSampleWorkspace(
        "Event",
        NumBanks=1,
        BankPixelWidth=1,
        NumEvents=0,
        SourceDistanceFromSample=14.1858856536088,
        BankDistanceFromSample=1.3,
    )
    clean_workspace(w)
    s = w.getSpectrum(0)
    for tof in starting_tof:
        s.addEventQuickly(float(tof), DateAndTime(0))

    # run correct_emission_time on workspace
    correct_frame.correct_emission_time(w)

    # convert the final tofs to wavelength and compare to expected values
    h = 6.62606896e-34
    m = 1.674927211e-27
    z = 15.4858856536088
    assert_allclose(w.getSpectrum(0).getTofs() * 10000 * h / (z * m), expected_wl, rtol=1e-4)


def test_correct_emission_time_30Hz(clean_workspace):
    # excepted wavelengths
    expected_wl = [
        2.55,
        2.65,
        2.75,
        2.85,
        2.95,
        3.05,
        3.15,
        3.25,
        3.35,
        3.45,
        3.55,
        3.65,
        3.75,
        3.85,
        3.95,
        4.05,
        4.15,
        4.25,
        4.35,
        4.45,
        4.55,
        4.65,
        4.75,
        4.85,
        4.95,
        5.05,
        5.15,
        5.25,
        5.35,
        5.45,
        5.55,
        5.65,
        5.75,
        5.85,
        5.95,
        6.05,
        6.15,
        9.75,
        9.85,
        9.95,
        10.05,
        10.15,
        10.25,
        10.35,
        10.45,
        10.55,
        10.65,
        10.75,
        10.85,
        10.95,
        11.05,
        11.15,
        11.25,
        11.35,
        11.45,
        11.55,
        11.65,
        11.75,
        11.85,
        11.95,
        12.05,
        12.15,
        12.25,
        12.35,
        12.45,
        12.55,
        12.65,
        12.75,
        12.85,
        12.95,
        13.05,
        13.15,
        13.25,
        13.35,
        13.45,
    ]
    starting_tof = [
        11814.9960422188,
        12273.7212567657,
        12732.4424149125,
        13191.1595610593,
        13649.8727396061,
        14108.5819949529,
        14567.2873714998,
        15025.9889136466,
        15484.6866657934,
        15943.3806723402,
        16402.070977687,
        16860.7576262338,
        17319.4406623807,
        17778.1201305275,
        18236.7960750743,
        18695.4685404211,
        19154.1375709679,
        19612.8032111147,
        20071.4655052616,
        20530.1244978084,
        20988.7802331552,
        21447.432755702,
        21906.0821098488,
        22364.7283399956,
        22823.3714905425,
        23282.0116058893,
        23740.6487304361,
        24199.2829085829,
        24657.9141847297,
        25116.5426032765,
        25575.1682086234,
        26033.7910451702,
        26492.411157317,
        26951.0285894638,
        27409.6433860106,
        27868.2555913574,
        28326.8652499043,
        44835.4913471896,
        45294.0379873364,
        45752.5837678833,
        46211.1287332301,
        46669.6729277769,
        47128.2163959237,
        47586.7591820705,
        48045.3013306173,
        48503.8428859642,
        48962.383892511,
        49420.9243946578,
        49879.4644368046,
        50338.0040633514,
        50796.5433186982,
        51255.0822472451,
        51713.6208933919,
        52172.1593015387,
        52630.6975160855,
        53089.2355814323,
        53547.7735419791,
        54006.3114421259,
        54464.8493262728,
        54923.3872388196,
        55381.9252241664,
        55840.4633267132,
        56299.00159086,
        56757.5400610068,
        57216.0787815537,
        57674.6177969005,
        58133.1571514473,
        58591.6968895941,
        59050.2370557409,
        59508.7776942877,
        59967.3188496346,
        60425.8605661814,
        60884.4028883282,
        61342.945860475,
        61801.4895270218,
    ]

    # Make a simple workspace with correct distances and add tofs to it
    w = CreateSampleWorkspace(
        "Event",
        NumBanks=1,
        BankPixelWidth=1,
        NumEvents=0,
        SourceDistanceFromSample=14.1395946855299,
        BankDistanceFromSample=4.0,
    )
    clean_workspace(w)
    s = w.getSpectrum(0)
    for tof in starting_tof:
        s.addEventQuickly(float(tof), DateAndTime(0))

    # run correct_emission_time on workspace
    correct_frame.correct_emission_time(w)

    # convert the final tofs to wavelength and compare to expected values
    h = 6.62606896e-34
    m = 1.674927211e-27
    z = 18.1395946855299
    assert_allclose(w.getSpectrum(0).getTofs() * 10000 * h / (z * m), expected_wl, rtol=1e-4)


@pytest.mark.parametrize(
    "wavelength, expected",
    [
        # λ < 2 Å branch
        (1.0, 68.555),
        # λ >= 2 Å branch
        (5.0, 126.578),
        # small positive λ — polynomial branch must not return a negative delay
        (0.1, "non_negative"),
        # zero wavelength returns 0 µs (linear segment)
        (0.0, 0.0),
        # negative wavelengths must raise ValueError
        (-1.0, ValueError),
        (-0.001, ValueError),
    ],
)
def test_emission_delay(wavelength, expected):
    if expected is ValueError:
        with pytest.raises(ValueError, match="wavelength must be non-negative"):
            correct_frame.emission_delay(wavelength)
    elif expected == "non_negative":
        assert correct_frame.emission_delay(wavelength) >= 0
    else:
        assert correct_frame.emission_delay(wavelength) == pytest.approx(expected, abs=1e-3)


if __name__ == "__main__":
    pytest.main([__file__])
