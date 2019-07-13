import pytest
from pytest import approx
from os.path import join as pj
from mantid.simpleapi import (Integration, SumSpectra)

from ornl.sans.sns.eqsans.normalisation import \
    (load_beam_flux_file, normalise_by_proton_charge_and_flux,
     normalise_by_flux)
from ornl.settings import amend_config, unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength)
from ornl.sans.samplelogs import SampleLogs


@pytest.fixture(scope='module')
def flux_file(refd):
    return pj(refd.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt')


@pytest.fixture(scope='module')
def data_ws(refd):
    with amend_config(data_dir=refd.new.eqsans):
        w = load_events('EQSANS_92353', output_workspace=uwd())
    return transform_to_wavelength(w, output_workspace=w.name())


def test_load_beam_flux_file(flux_file, data_ws):
    # No reference workspace
    w = load_beam_flux_file(flux_file)
    assert w.name().startswith('__')  # gets random hidden name
    i = Integration(w)
    assert i.dataY(0)[0] == approx(1)  # normalized
    assert max(w.dataY(0)) == approx(0.337, abs=0.001)
    # Reference workspace
    w = load_beam_flux_file(flux_file, ws_reference=data_ws)
    i = Integration(w)
    assert i.dataY(0)[0] == approx(0.79, abs=0.01)  # normalized
    assert max(w.dataY(0)) == approx(0.337, abs=0.001)
    assert w.dataX(0) == approx(data_ws.dataX(0))


def test_normalize_by_proton_charge_and_flux(flux_file, data_ws):
    flux_ws = load_beam_flux_file(flux_file, ws_reference=data_ws)
    w = normalise_by_proton_charge_and_flux(data_ws, flux_ws,
                                            output_workspace=uwd())
    pc = SampleLogs(data_ws).getProtonCharge()
    u = SumSpectra(w, OutputWorkspace=uwd()).dataY(0)
    u2 = SumSpectra(data_ws, OutputWorkspace=uwd()).dataY(0)
    assert u == approx(u2 / (flux_ws.dataY(0) * pc), rel=0.01)


def test_normalise_by_flux(flux_file, data_ws):
    w = normalise_by_flux(data_ws, flux_file, output_workspace='data_normed')
    assert w.name() == 'data_normed'
    u = SumSpectra(w, OutputWorkspace=uwd()).dataY(0)
    u2 = SumSpectra(data_ws, OutputWorkspace=uwd()).dataY(0)
    flux_ws = load_beam_flux_file(flux_file, ws_reference=data_ws)
    pc = SampleLogs(data_ws).getProtonCharge()
    assert u == approx(u2 / (flux_ws.dataY(0) * pc), rel=0.01)


if __name__ == '__main__':
    pytest.main()
