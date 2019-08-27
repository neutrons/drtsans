import pytest
import numpy as np

from ornl.sans.samplelogs import SampleLogs
from ornl.sans.hfir.normalisation import time, monitor

x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
z = (x + 5 * y) * 10 + 100


@pytest.mark.parametrize('generic_workspace',
                         [{'axis_values': x.tolist(),
                          'intensities': z.tolist()}],
                         indirect=True)
def test_normalization_by_time(generic_workspace):
    ws = generic_workspace
    print(ws.getNumberHistograms(), len(ws.dataX(0)), len(ws.dataY(0)))
    I_sam = z  # choose sample data from g or z
    I_sam_err = np.sqrt(z)

    t_sam = 5  # any value
    t_sam_err = 0.2  # 2 percent error

    I_samnorm = I_sam / t_sam
    I_norm_err = np.sqrt((I_sam_err / t_sam)**2 + (I_sam * t_sam_err / t_sam**2)**2)


def test_normalization_by_monitor():
    x, y = np.meshgrid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
    z = (x*0 + y*0 + 1)  # constant image
    sigma, mu = 0.3, 0  # gaussian image
    d = np.sqrt(x*x + y*y)
    g = np.exp(-((d-mu)**2/(2 * sigma**2)))
    z = (x + 5 * y) * 10 + 100

    I_sam = z  # choose sample data from g or z
    
    factor_is = 10**8
    flux_sam = 5 * 10**8
    I_samnorm = factor_is / flux_sam * I_sam


if __name__ == '__main__':
    pytest.main()