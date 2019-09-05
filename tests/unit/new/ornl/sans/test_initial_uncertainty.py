import numpy as np
import pytest
# https://docs.mantidproject.org/nightly/algorithms/LoadEmptyInstrument-v1.html
from mantid.simpleapi import LoadEmptyInstrument, Rebin, AddTimeSeriesLog, AddSampleLog
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
    wave_length_range = np.array([2.5, 6.5])  # A
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
        ws.dataX(i)[:] = wave_length_range  # A
        ws.dataY(i)[0] = intensity[i]
        ws.dataE(i)[0] = init_delta_intensity[i]
    # #### ABOVE THIS POINT WILL BE A TEST FIXTURE

    # Set uncertainties
    ws = set_init_uncertainties(ws, mask_band_gap=False)

    print('[TEST INFO] Workspace {} has {} spectra'.format(ws, ws.getNumberHistograms()))
    for ws_index in range(4):
        if np.isnan(gold_delta_intensity[ws_index]):
            assert np.isnan(ws.dataE(ws_index)[0])
        else:
            assert abs(ws.readE(ws_index)[0] -
                       gold_delta_intensity[ws_index]) < 1.E-10
        # END-IF
    # END-FOR

    # reuse the workspace to test with 'band gap'
    band_gap_test(ws)

    return


def band_gap_test(ws):
    """
    Test for the case with band gap
    :param ws: Workspace in Wavelength unit
    :return:
    """
    # Rebin to 2.5, 0.1, 6.5 and set value to intensity
    ws = Rebin(InputWorkspace=ws, Params='2.5, 0.1, 6.5')
    num_bins = ws.readY(0).shape[0]

    for iws in range(ws.getNumberHistograms()):
        vec_y_i = abs(np.random.randn(num_bins)) * 100
        vec_y_i[1] = 0
        vec_y_i[2] = np.nan
        ws.dataY(iws)[:] = vec_y_i[:]
    # END-FOR

    # add sample logs: chopper, lead_max, and skip_min
    AddTimeSeriesLog(ws, Name="BL6:Chop:Skf1:SpeedUserReq", Time="2010-01-01T00:00:00", Value=30.)
    AddTimeSeriesLog(ws, Name="BL6:Chop:Skf1:SpeedUserReq", Time="2010-01-01T00:30:00", Value=30.)

    AddSampleLog(Workspace=ws, LogName='lead_max', LogText='3.5', LogType='Number')
    AddSampleLog(Workspace=ws, LogName='skip_min', LogText='5.5', LogType='Number')

    # Set uncertainties
    ws = set_init_uncertainties(ws, mask_band_gap=True)

    # Check bins with wave length inside gap shall have uncertainties ZERO
    for iws in range(ws.getNumberHistograms()):
        assert ws.readE(iws)[1] == 1.
        assert np.isnan(ws.readE(iws)[2])
        assert np.allclose(ws.readE(iws)[10:30], np.zeros(20))
    # END-FOR

    return


if __name__ == '__main__':
    pytest.main()
