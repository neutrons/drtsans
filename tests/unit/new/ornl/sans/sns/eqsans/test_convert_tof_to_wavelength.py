import numpy as np
import pytest
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans.correct_frame import convert_to_wavelength


'''
# TODO expand on the next statement
# in master document section 3.3
# dev - Pete Peterson <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>
@pytest.mark.parametrize('generate_sans_generic_IDF',
                         [{'Nx': 2, 'Ny': 2}],
                         indirect=True)
@pytest.mark.parametrize('I, dI',
                         [(187., np.sqrt(187.)),
                          (0., 1.),
                          (np.nan, np.nan)],
                         ids=('I=187', 'I=0', 'I=NaN'))
def test_william(generate_sans_generic_IDF, I, dI):
    # generate a generic SANS instrument with a pixel of
    # the size and position specified in
    # sans-backend/documents/Master_document_022219.pdf
    with open(r'/tmp/GenericSANS_Definition.xml', 'w') as tmp:
        tmp.write(generate_sans_generic_IDF)
        tmp.close()
    ws = LoadEmptyInstrument(Filename=tmp.name, InstrumentName='GenericSANS')

    ws.dataY(1)[0] = I
    # TODO set uncertainty using non-existant function

    assert ws.dataY(1)[0] == I
    assert ws.dataE(1)[0] == dI
'''


TOF = [12345., 12346.]


# TODO expand on the next statement
# in master document section 3.3
# dev - Pete Peterson <petersonpf@ornl.gov>
# SME - Shuo Qian
@pytest.mark.parametrize('generic_workspace',
                         [{'dx': .005, 'dy': .004,
                           'zc': 2.5, 'L1': 10.1,
                           'axis_units': 'tof',
                           'axis_values': TOF,
                           'intensities': [[1., 1.], [1., 1.]]}],
                         indirect=True)
def test_shuo(generic_workspace):
    ws = generic_workspace  # friendly name
    samplelog = SampleLogs(ws)
    samplelog.insert('is_frame_skipping', False)

    ws = convert_to_wavelength(input_workspace=generic_workspace)

    specInfo = ws.spectrumInfo()
    source_sample = specInfo.l1()  # in meters

    # make sure the unit is wavelength
    assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
    # verify the individual wavelength values
    for i in range(4):
        sample_detector = specInfo.l2(i)  # to detector pixel in meters
        # equation taken
        lambda_exp = 3.9560346e-3 * np.array(TOF) / (source_sample + sample_detector)
        assert ws.dataX(i)[0] == pytest.approx(lambda_exp[0])
        assert ws.dataX(i)[1] == pytest.approx(lambda_exp[1])
        assert ws.dataX(i)[0] == pytest.approx(3.8760)


if __name__ == '__main__':
    pytest.main()
