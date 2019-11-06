import pytest
from pytest import approx
from os.path import join as pj

r"""
Hyperlinks to Mantid algorithms
SumSpectra <https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html>
"""
from mantid.simpleapi import SumSpectra

r"""
Hyperlinks to drtsans functions
amend_config, unique_workspace_dundername <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/settings.py>
SampleLogs <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/samplelogs.py>
load_events <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/load.py>
prepare_monitors <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/api.py>
normalize_by_time,...load_flux_to_monitor_ratio_file <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/tof/eqsans/nomalization.py>
"""  # noqa: E501
from drtsans.settings import amend_config, unique_workspace_dundername
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (load_events, transform_to_wavelength, prepare_monitors, normalize_by_flux,
                                normalize_by_proton_charge_and_flux, normalize_by_time, normalize_by_monitor)
from drtsans.tof.eqsans.normalization import load_beam_flux_file, load_flux_to_monitor_ratio_file

@pytest.fixture(scope='module')
def beam_flux(reference_dir):
    r"""Filepath to the flux file"""
    return pj(reference_dir.new.eqsans, 'test_normalization', 'beam_profile_flux.txt')

@pytest.fixture(scope='module')
def flux_to_monitor(reference_dir):
    r"""Filepath to the flux-to-monitor-ratio file"""
    return pj(reference_dir.new.eqsans, 'test_normalization', 'flux_to_monitor_ratio.nxs')


@pytest.fixture(scope='module')
def data_ws(reference_dir):
    r"""Two Mantid workspaces containing intensities versus wavelength for each of the EQSANS pixel-detectors.
    The two workspaces correspond to runs 92353 and 88565."""
    ws = dict()
    with amend_config(data_dir=reference_dir.new.eqsans):
        for run in ('92353', '88565'):
            w = load_events('EQSANS_{}'.format(run), output_workspace=unique_workspace_dundername())
            ws[run] = transform_to_wavelength(w, output_workspace=w.name())
    return ws


@pytest.fixture(scope='module')
def monitor_ws(reference_dir):
    r"""Single-spectrum Workspace containing wavelength-dependent monitor counts for run 88565."""
    ws = dict()
    with amend_config(data_dir=reference_dir.new.eqsans):
        for run in ('88565',):
            ws[run] = prepare_monitors(run)
    return ws


def test_load_beam_flux_file(beam_flux, data_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    Loads flux file
    /SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/sns/eqsans/test_normalization/beam_profile_flux.txt
    into a Mantid workspace.
    """
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_ws['92353'])
    assert flux_workspace.readY(0)[0] == approx(959270., abs=1.)
    assert max(flux_workspace.readY(0)) == approx(966276., abs=1.)
    assert flux_workspace.dataX(0) == approx(data_ws['92353'].dataX(0))


def test_normalize_by_proton_charge_and_flux(beam_flux, data_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)

    """
    data_workspace = data_ws['92353']  # intensities versus wavelength for run 92353
    # Load flux file
    # /SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/sns/eqsans/test_normalization/beam_profile_flux.txt
    # into a workspace
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_workspace)

    # Use drtsans normalizing function
    normalized_data_workspace = normalize_by_proton_charge_and_flux(data_workspace, flux_workspace,
                                                                    output_workspace=unique_workspace_dundername())

    # We run a simplified comparison. We merge all spectra of the individual pixel-detectors onto a single spectrum
    normalized_total_intensities = SumSpectra(normalized_data_workspace,
                                              OutputWorkspace=unique_workspace_dundername()).dataY(0)
    unnormalized_total_intensities = SumSpectra(data_workspace, OutputWorkspace=unique_workspace_dundername()).dataY(0)

    # Manually normalize the unnormalized_total_intensities and compare to the result from using drtsans
    # normalizing function
    good_proton_charge = SampleLogs(data_workspace).getProtonCharge()
    manual_normalized_intensities = unnormalized_total_intensities / (flux_workspace.readY(0) * good_proton_charge)

    # compare the two spectra don't deviate more than 1%.
    assert normalized_total_intensities == approx(manual_normalized_intensities, rel=0.01)


def test_load_flux_to_monitor_ratio_file(flux_to_monitor, data_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    Loads flux-to-monitor-ratio file
    /SNS/EQSANS/shared/sans-backend/data/new/ornl/sans/sns/eqsans/test_normalization/flux_to_monitor_ratio.nxs
    onto a mantid workspace
    """
    # Passing just the file to function load_flux_to_monitor_ratio_file will result in a workspace with the
    # wavelength binning as in the file.
    flux_to_monitor_workspace = load_flux_to_monitor_ratio_file(flux_to_monitor)
    # check that the workspace is a histogram (the number of wavelength boundaries is the number of ratios plus one)
    assert len(flux_to_monitor_workspace.dataX(0)) == 1 + len(flux_to_monitor_workspace.dataY(0))
    # check the number of wavelength bin boundaries is that of the input file.
    assert len(flux_to_monitor_workspace.dataX(0)) == 48664

    # Passing the file and a reference workspace to function load_flux_to_monitor_ratio_file will result
    # in a workspace with the wavelength binning as in the reference workspace.
    data_workspace = data_ws['88565']  # our reference workspace
    flux_to_monitor_workspace = load_flux_to_monitor_ratio_file(flux_to_monitor, data_workspace=data_workspace)
    # Check the wavelength bin boundaries are those of the reference workspace.
    assert flux_to_monitor_workspace.dataX(0) == approx(data_workspace.dataX(0), abs=1e-3)
    assert max(flux_to_monitor_workspace.dataY(0)) == approx(0.569, abs=1e-3)  # a simple check


def test_normalize_by_monitor(flux_to_monitor, data_ws, monitor_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # First we try normalization in frame-skipping mode, which should raise an exception
    data_workspace, monitor_workspace = data_ws['92353'], monitor_ws['88565']
    with pytest.raises(ValueError, match='not possible in frame-skipping'):
        # the below statement will raise an exception of class ValueError with an error message that
        # should contain the above "match" string
        data_workspace_normalized = normalize_by_monitor(data_workspace, flux_to_monitor, monitor_workspace,
                                                         output_workspace=unique_workspace_dundername())
    # Second we try normalization if non-skipping mode
    data_workspace, monitor_workspace = data_ws['88565'], monitor_ws['88565']  # data and monitor for run 88565
    # Use drtsans function to normalize by monitor counts
    data_workspace_normalized = normalize_by_monitor(data_workspace, flux_to_monitor, monitor_workspace,
                                                     output_workspace=unique_workspace_dundername())

    # Simplified test by comparing the intensity integrated over all pixel detectors and over all wavelength bins
    # after normalization
    # First we add the spectrum of all pixels into a single pixel
    data_workspace_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=data_workspace_normalized.name())
    # Second we integrate over all wavelength bins
    assert sum(data_workspace_normalized.dataY(0)) == approx(0.594, abs=1e-03)

    data_workspace_normalized.delete()  # clean-up from memory of temporary workspaces


def test_normalize_by_time(data_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    data_workspace = data_ws['92353']
    y, e = data_workspace.readY(42)[5], data_workspace.readE(42)[5]  # some meaningful choice

    normalized_data_workspace = normalize_by_time(data_workspace, output_workspace=unique_workspace_dundername())
    run_duration = SampleLogs(normalized_data_workspace).duration.value
    assert (y/run_duration, e/run_duration) == approx((normalized_data_workspace.readY(42)[5],
                                                       normalized_data_workspace.readE(42)[5]), abs=1e-6)
    assert SampleLogs(normalized_data_workspace).normalizing_duration.value == 'duration'
    normalized_data_workspace.delete()


def test_normalize_by_flux(beam_flux, flux_to_monitor, data_ws, monitor_ws):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # Normalize by flux and proton charge
    data_workspace = data_ws['92353']
    data_workspace_normalized = normalize_by_flux(data_workspace, beam_flux,
                                                  output_workspace=unique_workspace_dundername())
    summed_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=unique_workspace_dundername())
    # Compare to "manual" normalization
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_workspace)
    pc = SampleLogs(data_workspace).getProtonCharge()
    summed = SumSpectra(data_workspace, OutputWorkspace=unique_workspace_dundername())
    assert summed_normalized.readY(0) == approx(summed.readY(0) / (flux_workspace.readY(0) * pc), rel=0.01)
    [ws.delete() for ws in [data_workspace_normalized, flux_workspace, summed, summed_normalized]]

    # Normalize by monitor and flux-to-monitor ratio
    data_workspace, monitor_workspace = data_ws['88565'], monitor_ws['88565']
    data_workspace_normalized = normalize_by_flux(data_workspace, flux_to_monitor, method='monitor',
                                                  monitor_workspace=monitor_workspace,
                                                  output_workspace=unique_workspace_dundername())
    summed_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=unique_workspace_dundername())
    assert sum(summed_normalized.readY(0)) == approx(0.594, abs=1e-3)
    [ws.delete() for ws in [data_workspace_normalized, summed_normalized]]

    # Normalize by run duration
    data_workspace = data_ws['92353']
    duration = SampleLogs(data_workspace).duration.value
    total_intensity, total_intensity_error = sum(data_workspace.readY(42)), sum(data_workspace.readE(42))
    data_workspace_normalized = normalize_by_flux(data_workspace, 'duration', method='time',
                                                  output_workspace=unique_workspace_dundername())
    assert sum(data_workspace.readY(42)), sum(data_workspace.readE(42)) == approx((total_intensity / duration,
                                                                                   total_intensity_error / duration))
    data_workspace_normalized.delete()


if __name__ == '__main__':
    pytest.main([__file__])
