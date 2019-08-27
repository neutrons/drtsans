import numpy as np
import pytest
from ornl.sans.samplelogs import SampleLogs
from ornl.sans.sns.eqsans.correct_frame import convert_to_wavelength


def add_frame_skipping_log(ws):
    samplelog = SampleLogs(ws)
    samplelog.insert('is_frame_skipping', False)


# TODO expand on the next statement
# in master document section 3.3
# dev - Pete Peterson <petersonpf@ornl.gov>
# SME - William Heller <hellerwt@ornl.gov>
@pytest.mark.parametrize('generic_workspace',
                         [{'dx': .005, 'dy': .004,
                           'zc': 5., 'l1': 14.,
                           'axis_units': 'tof',
                           'axis_values': [15432.]}],
                         indirect=True)
def test_william(generic_workspace):
    ws = generic_workspace  # friendly name
    add_frame_skipping_log(ws)

    ws = convert_to_wavelength(input_workspace=generic_workspace)

    specInfo = ws.spectrumInfo()
    source_sample = specInfo.l1()  # in meters

    # make sure the unit is wavelength
    assert ws.getAxis(0).getUnit().caption() == 'Wavelength'
    # verify the individual wavelength values
    for i in range(4):
        sample_detector = specInfo.l2(i)  # to detector pixel in meters
        # equation taken
        lambda_exp = 3.9560346e-3 * np.array([15432.]) / (source_sample + sample_detector)
        assert ws.dataX(i)[0] == pytest.approx(lambda_exp[0])
        assert ws.dataX(i)[0] == pytest.approx(3.2131329446)


TOF = [12345., 12346.]


# TODO expand on the next statement
# in master document section 3.3
# dev - Pete Peterson <petersonpf@ornl.gov>
# SME - Shuo Qian
@pytest.mark.parametrize('generic_workspace',
                         [{'dx': .005, 'dy': .004,
                           'zc': 2.5, 'l1': 10.1,
                           'axis_units': 'tof',
                           'axis_values': TOF}],
                         indirect=True)
def test_shuo(generic_workspace):
    ws = generic_workspace  # friendly name
    add_frame_skipping_log(ws)

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
        assert ws.dataX(i)[0] == pytest.approx(3.875969)  # Shuo asked for 3.8760


if __name__ == '__main__':
    pytest.main()
