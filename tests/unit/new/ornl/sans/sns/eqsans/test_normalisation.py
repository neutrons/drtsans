import pytest
from pytest import approx
from os.path import join as pj
from mantid.simpleapi import (Integration, SumSpectra)

from ornl.sans.sns.eqsans.normalisation import \
    (load_beam_flux_file, normalise_by_proton_charge_and_flux,
     load_flux_to_monitor_ratio_file, normalise_by_monitor, normalise_by_flux)
from ornl.settings import amend_config, unique_workspace_dundername as uwd
from ornl.sans.sns.eqsans import (load_events, transform_to_wavelength,
                                  prepare_monitors)
from ornl.sans.samplelogs import SampleLogs


@pytest.fixture(scope='module')
def beam_flux(refd):
    return pj(refd.new.eqsans, 'test_normalisation', 'beam_profile_flux.txt')


@pytest.fixture(scope='module')
def flux_to_monitor(refd):
    return pj(refd.new.eqsans,
              'test_normalisation', 'flux_to_monitor_ratio.nxs')


@pytest.fixture(scope='module')
def data_ws(refd):
    ws = dict()
    with amend_config(data_dir=refd.new.eqsans):
        for run in ('92353', '88565'):
            w = load_events('EQSANS_{}'.format(run), output_workspace=uwd())
            ws[run] = transform_to_wavelength(w, output_workspace=w.name())
    return ws


@pytest.fixture(scope='module')
def monitor_ws(refd):
    ws = dict()
    with amend_config(data_dir=refd.new.eqsans):
        for run in ('88565',):
            ws[run] = prepare_monitors(run)
    return ws


def test_load_beam_flux_file(beam_flux, data_ws):
    # No reference workspace
    w = load_beam_flux_file(beam_flux)
    assert w.name().startswith('__')  # gets random hidden name
    i = Integration(w)
    assert i.dataY(0)[0] == approx(1)  # normalized
    assert max(w.dataY(0)) == approx(0.337, abs=0.001)
    # Reference workspace
    w = load_beam_flux_file(beam_flux, ws_reference=data_ws['92353'])
    i = Integration(w)
    assert i.dataY(0)[0] == approx(0.79, abs=0.01)  # normalized
    assert max(w.dataY(0)) == approx(0.337, abs=0.001)
    assert w.dataX(0) == approx(data_ws['92353'].dataX(0))


def test_normalize_by_proton_charge_and_flux(beam_flux, data_ws):
    dws = data_ws['92353']
    flux_ws = load_beam_flux_file(beam_flux, ws_reference=dws)
    w = normalise_by_proton_charge_and_flux(dws, flux_ws,
                                            output_workspace=uwd())
    pc = SampleLogs(dws).getProtonCharge()
    u = SumSpectra(w, OutputWorkspace=uwd()).dataY(0)
    u2 = SumSpectra(dws, OutputWorkspace=uwd()).dataY(0)
    assert u == approx(u2 / (flux_ws.dataY(0) * pc), rel=0.01)


def test_normalise_by_flux(beam_flux, data_ws):
    dws = data_ws['92353']
    w = normalise_by_flux(dws, beam_flux, output_workspace='data_normed')
    assert w.name() == 'data_normed'
    u = SumSpectra(w, OutputWorkspace=uwd()).dataY(0)
    u2 = SumSpectra(dws, OutputWorkspace=uwd()).dataY(0)
    flux_ws = load_beam_flux_file(beam_flux, ws_reference=dws)
    pc = SampleLogs(dws).getProtonCharge()
    assert u == approx(u2 / (flux_ws.dataY(0) * pc), rel=0.01)


def test_load_flux_to_monitor_ratio_file(flux_to_monitor, data_ws):
    # No reference workspace
    w = load_flux_to_monitor_ratio_file(flux_to_monitor)
    assert len(w.dataX(0)) == 1 + len(w.dataY(0))
    assert len(w.dataX(0)) == 3283
    # Reference workspace
    dws = data_ws['88565']
    w = load_flux_to_monitor_ratio_file(flux_to_monitor, ws_reference=dws)
    assert w.dataX(0) == approx(dws.dataX(0), abs=1e-3)
    assert max(w.dataY(0)) * 1e7 == approx(4.25, abs=0.01)


def test_normalise_by_monitor(flux_to_monitor, data_ws, monitor_ws):
    dws, mon = data_ws['88565'], monitor_ws['88565']
    w = normalise_by_monitor(dws, mon, flux_to_monitor,
                             output_workspace=uwd())
    w.delete()


if __name__ == '__main__':
    pytest.main()
