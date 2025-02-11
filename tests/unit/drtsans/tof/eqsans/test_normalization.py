import pytest
from pytest import approx
from os.path import join as pj

r"""
Hyperlinks to Mantid algorithms
SumSpectra <https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html>
amend_config <https://docs.mantidproject.org/nightly/api/python/mantid/kernel/AmendConfig.html>
"""
from mantid.simpleapi import SumSpectra, mtd, Scale, CloneWorkspace
from mantid.kernel import amend_config

r"""
Hyperlinks to drtsans functions
SampleLogs <https://github.com/neutrons/drtsans/blob/next/src/drtsans/samplelogs.py>
load_events <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/load.py>
prepare_monitors <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/api.py>
normalize_by_time,...load_flux_to_monitor_ratio_file <https://github.com/neutrons/drtsans/blob/next/src/drtsans/tof/eqsans/nomalization.py>
"""  # noqa: E501
from drtsans.samplelogs import SampleLogs
from drtsans.tof.eqsans import (
    load_events,
    transform_to_wavelength,
    set_init_uncertainties,
    prepare_monitors,
    normalize_by_flux,
    normalize_by_proton_charge_and_flux,
    normalize_by_time,
    normalize_by_monitor,
    ZeroNeutronFluxError,
)
from drtsans.tof.eqsans.normalization import (
    load_beam_flux_file,
    load_flux_to_monitor_ratio_file,
)


@pytest.fixture(scope="module")
def beam_flux(datarepo_dir):
    r"""Filepath to the flux file"""
    return pj(datarepo_dir.eqsans, "test_normalization", "beam_profile_flux.txt")


@pytest.fixture(scope="module")
def flux_to_monitor(datarepo_dir):
    r"""Filepath to the flux-to-monitor-ratio file"""
    return pj(datarepo_dir.eqsans, "test_normalization", "flux_to_monitor_ratio.nxs")


@pytest.fixture(scope="module")
def data_ws(datarepo_dir):
    r"""Two Mantid workspaces containing intensities versus wavelength for each of the EQSANS pixel-detectors.
    The two workspaces correspond to runs 92353 and 88565."""
    ws = dict()
    with amend_config(data_dir=datarepo_dir.eqsans):
        for run in ("92353", "88565"):
            w = load_events(
                f"EQSANS_{run}.nxs.h5",
                output_workspace=mtd.unique_hidden_name(),
            )
            ws[run], bands = transform_to_wavelength(w, output_workspace=w.name())
            ws[run] = set_init_uncertainties(ws[run])
    return ws


@pytest.fixture(scope="module")
def monitor_ws(datarepo_dir):
    r"""Single-spectrum Workspace containing wavelength-dependent monitor counts for run 88565."""
    ws = dict()
    with amend_config(facility="SNS", instrument="EQSANS", data_dir=datarepo_dir.eqsans):
        for run in ("88565",):
            ws[run] = prepare_monitors(run)
    return ws


@pytest.mark.datarepo
def test_load_beam_flux_file(beam_flux, data_ws, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    Loads flux file
    /SNS/EQSANS/shared/sans-backend/data/ornl/sans/sns/eqsans/test_normalization/beam_profile_flux.txt
    into a Mantid workspace.
    """
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_ws["92353"])
    clean_workspace(flux_workspace)
    assert flux_workspace.readY(0)[0] == approx(959270.0, abs=1.0)
    assert max(flux_workspace.readY(0)) == approx(966276.0, abs=1.0)
    assert flux_workspace.dataX(0) == approx(data_ws["92353"].dataX(0))


@pytest.mark.datarepo
def test_normalize_by_proton_charge_and_flux(beam_flux, data_ws, temp_workspace_name, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)

    """
    data_workspace = data_ws["92353"]  # intensities versus wavelength for run 92353
    # Load flux file
    # /SNS/EQSANS/shared/sans-backend/data/ornl/sans/sns/eqsans/test_normalization/beam_profile_flux.txt
    # into a workspace
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_workspace)
    clean_workspace(flux_workspace)

    # Use drtsans normalizing function
    normalized_data_workspace = normalize_by_proton_charge_and_flux(
        data_workspace, flux_workspace, output_workspace=temp_workspace_name()
    )

    # We run a simplified comparison. We merge all spectra of the individual pixel-detectors onto a single spectrum
    normalized_total_intensities = SumSpectra(normalized_data_workspace, OutputWorkspace=temp_workspace_name()).dataY(
        0
    )
    unnormalized_total_intensities = SumSpectra(data_workspace, OutputWorkspace=temp_workspace_name()).dataY(0)

    # Manually normalize the unnormalized_total_intensities and compare to the result from using drtsans
    # normalizing function
    good_proton_charge = SampleLogs(data_workspace).getProtonCharge()
    manual_normalized_intensities = unnormalized_total_intensities / (flux_workspace.readY(0) * good_proton_charge)

    # compare the two spectra don't deviate more than 1%.
    assert normalized_total_intensities == approx(manual_normalized_intensities, rel=0.01)


@pytest.mark.datarepo
def test_normalize_by_proton_charge_and_flux_no_proton_charge(
    beam_flux, data_ws, temp_workspace_name, clean_workspace
):
    r"""
    (This test was introduced prior to the testset with the instrument team)

    """
    data_workspace = CloneWorkspace(data_ws["92353"])  # intensities versus wavelength for run 92353
    SampleLogs(data_workspace).insert("gd_prtn_chrg", 0.0)  # set proton charge to zero
    # Load flux file
    # /SNS/EQSANS/shared/sans-backend/data/ornl/sans/sns/eqsans/test_normalization/beam_profile_flux.txt
    # into a workspace
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_workspace)
    clean_workspace(flux_workspace)
    output_workspace = temp_workspace_name()

    # Use drtsans normalizing function
    with pytest.raises(ZeroNeutronFluxError) as except_info:
        normalize_by_proton_charge_and_flux(data_workspace, flux_workspace, output_workspace)

    assert f"Zero neutron flux for workspace: {output_workspace}" in str(except_info.value)


@pytest.mark.datarepo
def test_load_flux_to_monitor_ratio_file(flux_to_monitor, data_ws, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    Loads flux-to-monitor-ratio file
    /SNS/EQSANS/shared/sans-backend/data/ornl/sans/sns/eqsans/test_normalization/flux_to_monitor_ratio.nxs
    onto a mantid workspace
    """
    # Passing just the file to function load_flux_to_monitor_ratio_file will result in a workspace with the
    # wavelength binning as in the file.
    flux_to_monitor_workspace = load_flux_to_monitor_ratio_file(flux_to_monitor)
    clean_workspace(flux_to_monitor_workspace)
    # check that the workspace is a histogram (the number of wavelength boundaries is the number of ratios plus one)
    assert len(flux_to_monitor_workspace.dataX(0)) == 1 + len(flux_to_monitor_workspace.dataY(0))
    # check the number of wavelength bin boundaries is that of the input file.
    assert len(flux_to_monitor_workspace.dataX(0)) == 48664

    # Passing the file and a reference workspace to function load_flux_to_monitor_ratio_file will result
    # in a workspace with the wavelength binning as in the reference workspace.
    data_workspace = data_ws["88565"]  # our reference workspace
    flux_to_monitor_workspace = load_flux_to_monitor_ratio_file(flux_to_monitor, data_workspace=data_workspace)
    clean_workspace(flux_to_monitor_workspace)
    # Check the wavelength bin boundaries are those of the reference workspace.
    assert flux_to_monitor_workspace.dataX(0) == approx(data_workspace.dataX(0), abs=1e-3)
    assert max(flux_to_monitor_workspace.dataY(0)) == approx(0.569, abs=1e-3)  # a simple check


@pytest.mark.datarepo
def test_normalize_by_monitor(flux_to_monitor, data_ws, monitor_ws, temp_workspace_name):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # First we try normalization in frame-skipping mode, which should raise an exception
    data_workspace, monitor_workspace = data_ws["92353"], monitor_ws["88565"]
    with pytest.raises(ValueError, match="not possible in frame-skipping"):
        # the below statement will raise an exception of class ValueError with an error message that
        # should contain the above "match" string
        data_workspace_normalized = normalize_by_monitor(
            data_workspace,
            flux_to_monitor,
            monitor_workspace,
            output_workspace=temp_workspace_name(),
        )
    # Second we try normalization if non-skipping mode
    data_workspace, monitor_workspace = (
        data_ws["88565"],
        monitor_ws["88565"],
    )  # data and monitor for run 88565
    # Use drtsans function to normalize by monitor counts
    data_workspace_normalized = normalize_by_monitor(
        data_workspace,
        flux_to_monitor,
        monitor_workspace,
        output_workspace=temp_workspace_name(),
    )

    # Simplified test by checking the intensity integrated over all pixel detectors and over all wavelength bins
    # after normalization
    # First we add the spectrum of all pixels into a single pixel
    data_workspace_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=data_workspace_normalized.name())
    # Second we integrate over all wavelength bins and check the value  will not change as the code in the
    # repository evolves
    assert sum(data_workspace_normalized.dataY(0)) == approx(0.621, abs=1e-03)


@pytest.mark.datarepo
def test_normalize_by_monitor_zero_counts(flux_to_monitor, data_ws, monitor_ws, temp_workspace_name):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    # First we try normalization in frame-skipping mode, which should raise an exception
    data_workspace, monitor_workspace = data_ws["92353"], monitor_ws["88565"]
    with pytest.raises(ValueError, match="not possible in frame-skipping"):
        # the below statement will raise an exception of class ValueError with an error message that
        # should contain the above "match" string
        normalize_by_monitor(
            data_workspace,
            flux_to_monitor,
            monitor_workspace,
            output_workspace=temp_workspace_name(),
        )
    # Second we try normalization if non-skipping mode
    data_workspace = data_ws["88565"]
    monitor_workspace = CloneWorkspace(monitor_ws["88565"])
    monitor_workspace = Scale(monitor_workspace, 0.0)  # set monitor counts to zero
    output_workspace = temp_workspace_name()
    # data and monitor for run 88565
    # Use drtsans function to normalize by monitor counts
    with pytest.raises(ZeroNeutronFluxError) as except_info:
        normalize_by_monitor(
            data_workspace,
            flux_to_monitor,
            monitor_workspace,
            output_workspace,
        )

    assert f"Zero neutron flux for workspace: {output_workspace}" in str(except_info.value)


@pytest.mark.datarepo
def test_normalize_by_time(data_ws, temp_workspace_name):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    """
    data_workspace = data_ws["92353"]  # intensities versus wavelength for run 92353

    # use drtsans normalizing function
    data_workspace_normalized = normalize_by_time(data_workspace, output_workspace=temp_workspace_name())

    # check we looked log entry 'duration' in order to find out the time duration of the run
    assert SampleLogs(data_workspace_normalized).normalizing_duration.value == "duration"
    run_duration = SampleLogs(data_workspace_normalized).duration.value  # run duration, in seconds
    assert run_duration == pytest.approx(101.8, abs=0.1)

    # Simplified test by checking the intensity integrated over all pixel detectors and over all wavelength bins
    # after normalization
    # First we add the spectrum of all pixels into a single pixel
    data_workspace_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=data_workspace_normalized.name())
    # Second we integrate over all wavelength bins and check the value will not change as the code in the repository
    # evolves
    assert sum(data_workspace_normalized.dataY(0)) == approx(2560.5, abs=1.0)


@pytest.mark.datarepo
def test_normalize_by_flux(beam_flux, flux_to_monitor, data_ws, monitor_ws, temp_workspace_name, clean_workspace):
    r"""
    (This test was introduced prior to the testset with the instrument team)
    Function normalize_by_flux is a front to the three time of normalization we can carry out.
    """
    #
    # First we normalize by flux and proton charge with method='proton charge'
    #
    data_workspace = data_ws["92353"]
    data_workspace_normalized = normalize_by_flux(
        data_workspace,
        beam_flux,
        method="proton charge",
        output_workspace=temp_workspace_name(),
    )
    # we carry a simplified test whereby we will sum all pixel-detector spectra into a single spectrum
    summed_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=temp_workspace_name())
    summed_normalized_intensities = summed_normalized.readY(0)  # there's only one spectrum, that with index 0

    # Compare the output of calling function "normalize_by_flux" to a "manual" normalization by carrying out the
    # individual normalizing steps one after the other.
    flux_workspace = load_beam_flux_file(beam_flux, data_workspace=data_workspace)  # first load the flux file
    clean_workspace(flux_workspace)
    proton_charge = SampleLogs(data_workspace).getProtonCharge()  # find the proton charge
    summed = SumSpectra(data_workspace, OutputWorkspace=temp_workspace_name())
    manual_summed_normalized_intensities = summed.readY(0) / (flux_workspace.readY(0) * proton_charge)

    # compare now output of calling function "normalize_by_flux" to the "manual" normalization
    assert summed_normalized_intensities == pytest.approx(manual_summed_normalized_intensities, rel=0.001)

    #
    # Second we normalize by monitor and flux-to-monitor ratio with method='monitor'
    #
    data_workspace, monitor_workspace = data_ws["88565"], monitor_ws["88565"]
    data_workspace_normalized = normalize_by_flux(
        data_workspace,
        flux_to_monitor,
        method="monitor",
        monitor_workspace=monitor_workspace,
        output_workspace=temp_workspace_name(),
    )
    # we carry a simplified test whereby we will sum all pixel-detector spectra into a single spectrum
    summed_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=temp_workspace_name())
    # then we integrate this single spectrum over all wavelengths
    total_normalized_intensity = sum(summed_normalized.readY(0))
    # here we just check that the result will not change as the code in the repository evolves
    assert total_normalized_intensity == approx(0.621, abs=1e-3)

    #
    # Third we normalize by run duration with method='time'
    #
    data_workspace = data_ws["92353"]  # intensities versus wavelength for run 92353
    # call normalization by time using the log entry 'duration' when searching for the duration of the run
    data_workspace_normalized = normalize_by_flux(
        data_workspace,
        "duration",
        method="time",
        output_workspace=temp_workspace_name(),
    )

    # check we looked log entry 'duration' in order to find out the time duration of the run
    assert SampleLogs(data_workspace_normalized).normalizing_duration.value == "duration"
    run_duration = SampleLogs(data_workspace_normalized).duration.value  # run duration, in seconds
    assert run_duration == pytest.approx(101.8, abs=0.1)

    # Simplified test by checking the intensity integrated over all pixel detectors and over all wavelength bins
    # after normalization
    # First we add the spectrum of all pixels into a single pixel
    data_workspace_normalized = SumSpectra(data_workspace_normalized, OutputWorkspace=data_workspace_normalized.name())
    # Second we integrate over all wavelength bins and check the value will not change as the code in the repository
    # evolves
    assert sum(data_workspace_normalized.dataY(0)) == approx(2560, abs=1.0)


if __name__ == "__main__":
    pytest.main([__file__])
