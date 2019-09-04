import numpy as np
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
from mantid.simpleapi import LoadEmptyInstrument
from ornl.sans.process_uncertainties import set_init_uncertainties


# This implements Issue #163
# in master document section 3.3
# dev - Pete Peterson <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>
@pytest.mark.parametrize('generic_IDF',
                         [{'Nx': 2, 'Ny': 2,
                           'dx': 0.005, 'dy': 0.004, 'zc': 2.5}],
                         indirect=True)
def test_initial_uncertainty(generic_IDF):
    """
    Test initial uncertainty after histogram data is converted to unit
    wavelength for a TOF instrument
    :param generic_IDF,: IDF to generate
    :return:i
    """
    # Range of TOF
    TOF = np.array([12345., 12346.])  # microseconds
    intensity = np.array([0., 187., 1., np.nan])
    init_delta_intensity = np.random.randn(intensity.shape[0], )
    gold_delta_intensity = np.array([1., np.sqrt(187.), 1., np.nan])

    # Generate a generic SANS instrument with a pixel of
    # the size and position specified in
    # sans-backend/documents/Master_document_022219.pdf
    with open(r'/tmp/GenericSANS2_Definition.xml', 'w') as tmp:
        tmp.write(generic_IDF)
        tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS2',
                             OutputWorkspace='test_uncertainty')
    ws.getAxis(0).setUnit('Wavelength')
    # assume that the TOF is already frame corrected
    for i in range(4):
        ws.dataX(i)[:] = TOF  # microseconds
        ws.dataY(i)[0] = intensity[i]
        ws.dataE(i)[0] = init_delta_intensity[i]
    # #### ABOVE THIS POINT WILL BE A TEST FIXTURE

    # Set uncertainties
    ws = set_init_uncertainties(ws)

    print('[INFO] Workspace {} has {} spectra'.format(ws, ws.getNumberHistograms()))
    for ws_index in range(4):
        if np.isnan(gold_delta_intensity[ws_index]):
            assert np.isnan(ws.dataE(ws_index)[0])
        else:
            assert abs(ws.dataE(ws_index)[0] -
                       gold_delta_intensity[ws_index]) < 1.E-10
        # END-IF
    # END-FOR

    # TODO - reuse the workspace to test with 'band gap'
    band_gap_test(ws)

    return


def band_gap_test(ws):

    # 1. Rebin to 2.5, 0.1, 6.5

    # 2. add sample logs

    # 3. Set uncertainties

    # 4. Check bins with wave length inside gap shall have uncertainties ZERO

    return


def _set_eqsans_ws(ws):
    """
    Set properties to workspace including ...
    :param ws:
    :return:
    """
    # TODO - Implement!!!

if __name__ == '__main__':
    pytest.main()
