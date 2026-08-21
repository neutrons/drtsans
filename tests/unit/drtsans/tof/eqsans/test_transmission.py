import numpy as np
import pytest
from os.path import join as pjn
from numpy.testing import assert_almost_equal

# https://docs.mantidproject.org/nightly/algorithms/LoadNexus-v1.html
from mantid.simpleapi import CreateWorkspace, LoadNexus, mtd

# https://github.com/neutrons/drtsans/blob/next/src/drtsans/settings.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/correct_frame.py
# https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/transmission.py
from drtsans.settings import namedtuplefy  # noqa: E402
from drtsans.samplelogs import SampleLogs  # noqa: E402
from drtsans.tof.eqsans.correct_frame import TransmittedBands, transmitted_bands  # noqa: E402
from drtsans.tof.eqsans.transmission import fit_band, fit_raw_transmission  # noqa: E402
from drtsans.tof.eqsans.geometry import beam_radius  # noqa: E402
from drtsans.wavelength import Wband  # noqa: E402
import drtsans.tof.eqsans.transmission as transmission_module  # noqa: E402


@pytest.fixture(scope="module")
@namedtuplefy
def trasmission_data(datarepo_dir):
    data_dir = pjn(datarepo_dir.eqsans, "test_transmission")
    a = LoadNexus(pjn(data_dir, "raw_transmission.nxs"))
    b = LoadNexus(pjn(data_dir, "sample.nxs"))
    c = LoadNexus(pjn(data_dir, "raw_transmission_skip.nxs"))
    d = LoadNexus(pjn(data_dir, "sample_skip.nxs"))
    for workspace in (a, b, c, d):
        sample_logs = SampleLogs(workspace)
        sample_logs.insert("low_tof_clip", 0.0, unit="ms")
        sample_logs.insert("low_tof_clip", 0.0, unit="ms")
    return dict(data_dir=data_dir, raw=a, sample=b, raw_skip=c, sample_skip=d)


@pytest.mark.datarepo
def test_beam_radius(trasmission_data):
    r"""Verify the beam radius is correctly calculated using the source and sample apertures"""
    assert_almost_equal(beam_radius(trasmission_data.sample), 10.575980, decimal=6)
    assert_almost_equal(beam_radius(trasmission_data.sample_skip), 22.156863, decimal=6)


@pytest.mark.datarepo
def test_fit_band(trasmission_data, clean_workspace):
    r"""
    Verify the fitting of the raw transmissions over wavelength bands provides same fitting results.

    This test doesn't verify the correctness of the fit, only that the fit doesn't inadvertently change after
    changes to the source code are introduced. We use the goodness of fit (chi-square) to assess no changes.
    """
    # Non-skip mode
    bands = transmitted_bands(trasmission_data.raw)  # obtain the wavelength band
    ws_output, mantid_fit_output = fit_band(trasmission_data.raw, bands.lead)
    clean_workspace(ws_output)
    clean_workspace(mantid_fit_output.OutputWorkspace)
    clean_workspace(mantid_fit_output.OutputParameters)
    clean_workspace(mantid_fit_output.OutputNormalisedCovarianceMatrix)
    assert_almost_equal(mantid_fit_output.OutputChi2overDoF, 1.1, decimal=1)

    # Frame-skipping mode
    bands = transmitted_bands(trasmission_data.raw_skip)  # obtain the lead and skipped wavelength bands
    # fit raw transmission values in the wavelength band corresponding to the lead pulse
    ws_output, mantid_fit_output = fit_band(trasmission_data.raw_skip, bands.lead)
    clean_workspace(ws_output)
    clean_workspace(mantid_fit_output.OutputWorkspace)
    clean_workspace(mantid_fit_output.OutputParameters)
    clean_workspace(mantid_fit_output.OutputNormalisedCovarianceMatrix)
    assert_almost_equal(mantid_fit_output.OutputChi2overDoF, 1.1, decimal=1)
    # fit raw transmission values in the wavelength band corresponding to the skipped pulse
    ws_output, mantid_fit_output = fit_band(trasmission_data.raw_skip, bands.skip)
    clean_workspace(ws_output)
    clean_workspace(mantid_fit_output.OutputWorkspace)
    clean_workspace(mantid_fit_output.OutputParameters)
    clean_workspace(mantid_fit_output.OutputNormalisedCovarianceMatrix)
    assert_almost_equal(mantid_fit_output.OutputChi2overDoF, 3.6, decimal=0)


@pytest.mark.datarepo
def test_fit_raw(trasmission_data, clean_workspace):
    r"""
    Verify the fitting of raw transmission workspaces provides same fitting results.

    This test doesn't verify the correctness of the fit, only that the fit doesn't inadvertently change after
    changes to the source code are introduced. We use the goodness of fit (chi-square) to assess no changes.
    """
    # Non-skip mode
    fitting_results = fit_raw_transmission(trasmission_data.raw, output_workspace=mtd.unique_hidden_name())
    clean_workspace(fitting_results.transmission)
    clean_workspace(fitting_results.lead_transmission)
    clean_workspace(fitting_results.lead_mantid_fit.OutputWorkspace)
    clean_workspace(fitting_results.lead_mantid_fit.OutputNormalisedCovarianceMatrix)
    clean_workspace(fitting_results.lead_mantid_fit.OutputParameters)
    assert_almost_equal(fitting_results.lead_mantid_fit.OutputChi2overDoF, 1.1, decimal=1)

    # Frame-skipping mode
    fitting_results = fit_raw_transmission(trasmission_data.raw_skip, output_workspace=mtd.unique_hidden_name())
    clean_workspace(fitting_results.transmission)
    clean_workspace(fitting_results.lead_transmission)
    clean_workspace(fitting_results.lead_mantid_fit.OutputWorkspace)
    clean_workspace(fitting_results.lead_mantid_fit.OutputNormalisedCovarianceMatrix)
    clean_workspace(fitting_results.lead_mantid_fit.OutputParameters)
    clean_workspace(fitting_results.skip_transmission)
    clean_workspace(fitting_results.skip_mantid_fit.OutputWorkspace)
    clean_workspace(fitting_results.skip_mantid_fit.OutputNormalisedCovarianceMatrix)
    clean_workspace(fitting_results.skip_mantid_fit.OutputParameters)
    assert_almost_equal(fitting_results.lead_mantid_fit.OutputChi2overDoF, 1.1, decimal=1)
    assert_almost_equal(fitting_results.skip_mantid_fit.OutputChi2overDoF, 1.4, decimal=1)


def test_fit_band_bypasses_fit_for_single_valid_point(clean_workspace):
    input_workspace = CreateWorkspace(
        DataX=[2.45, 2.55],
        DataY=[0.873],
        DataE=[0.012],
        UnitX="Wavelength",
        OutputWorkspace=mtd.unique_hidden_name(),
    )

    output_workspace, mantid_fit_output = fit_band(
        input_workspace,
        Wband(2.45, 2.55),
        output_workspace=mtd.unique_hidden_name(),
    )

    assert mantid_fit_output is None
    assert output_workspace.readY(0).tolist() == pytest.approx([0.873])
    assert output_workspace.readE(0).tolist() == pytest.approx([0.012])
    clean_workspace(input_workspace)
    clean_workspace(output_workspace)


def test_fit_band_rejects_empty_band(clean_workspace):
    input_workspace = CreateWorkspace(
        DataX=[2.45, 2.55],
        DataY=[np.nan],
        DataE=[np.nan],
        UnitX="Wavelength",
        OutputWorkspace=mtd.unique_hidden_name(),
    )

    with pytest.raises(RuntimeError, match="No valid transmission points"):
        fit_band(input_workspace, Wband(2.45, 2.55), output_workspace=mtd.unique_hidden_name())
    clean_workspace(input_workspace)


def test_fit_raw_handles_lead_and_skip_bands_independently(monkeypatch, clean_workspace):
    input_workspace = CreateWorkspace(
        DataX=[1.0, 2.0, 3.0, 4.0],
        DataY=[0.5, 0.6, 0.7],
        DataE=[0.01, 0.01, 0.01],
        UnitX="Wavelength",
        OutputWorkspace=mtd.unique_hidden_name(),
    )
    monkeypatch.setattr(
        transmission_module,
        "transmitted_bands_clipped",
        lambda _: TransmittedBands(Wband(1.0, 2.0), Wband(2.0, 4.0)),
    )

    fitting_results = transmission_module.fit_raw_transmission(
        input_workspace,
        output_workspace=mtd.unique_hidden_name(),
    )

    assert fitting_results.lead_mantid_fit is None
    assert fitting_results.skip_mantid_fit is not None
    assert fitting_results.lead_transmission.readY(0).tolist() == pytest.approx([0.5, 0.0, 0.0])
    clean_workspace(input_workspace)
    clean_workspace(fitting_results.transmission)
    clean_workspace(fitting_results.lead_transmission)
    clean_workspace(fitting_results.skip_transmission)
    clean_workspace(fitting_results.skip_mantid_fit.OutputWorkspace)
    clean_workspace(fitting_results.skip_mantid_fit.OutputNormalisedCovarianceMatrix)
    clean_workspace(fitting_results.skip_mantid_fit.OutputParameters)


if __name__ == "__main__":
    pytest.main([__file__])
