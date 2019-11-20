from drtsans.plots import plot_IQmod, plot_IQazimuthal
from drtsans.dataobjects import IQmod, IQazimuthal
import numpy as np
import os
import pytest


def fileCheckAndRemove(filename, remove=True):
    assert os.path.exists(filename), '"{}" does not exist'.format(filename)
    assert os.path.getsize(filename) > 100, '"{}" is too small'.format(filename)
    if remove:
        os.remove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQmod.png'),
                          ('d3', 'test_IQmod.json')],
                         ids=['mpl', 'd3'])
def test_IQmod(backend, filename):
    x = np.linspace(0., 4*np.pi, 50)
    e = np.zeros(50) + .1
    data = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    plot_IQmod([data], filename=filename, backend=backend)

    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQmod_multi.png'),
                          ('d3', 'test_IQmod_multi.json')],
                         ids=['mpl', 'd3'])
def test_IQmod_multi(backend, filename):
    x = np.linspace(0., 4*np.pi, 50)
    e = np.zeros(50) + .1
    data1 = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    data2 = IQmod(intensity=np.cos(x), error=e, mod_q=x)

    plot_IQmod([data1, data2], filename=filename, backend=backend)
    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQazimuthal_1d.png'),
                          ('d3', 'test_IQazimuthal_1d.json')],
                         ids=['mpl', 'd3'])
def test_IQazimuthal_1d(backend, filename):
    x = np.linspace(0., 4*np.pi, 50)
    y = np.linspace(.5 * np.pi, 4.5*np.pi, 50)
    error = np.zeros((50, 50))
    data = error[:]
    for i in range(x.size):
        for j in range(y.size):
            data[i, j] = np.sin(x[i]) + np.cos(y[j])
    data = IQazimuthal(intensity=data, error=error, qx=x, qy=y)

    plot_IQazimuthal(data, filename=filename, backend=backend)
    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQazimuthal_2d.png'),
                          ('d3', 'test_IQazimuthal_2d.json')],
                         ids=['mpl', 'd3'])
def test_IQazimuthal_2d(backend, filename):
    x, y = np.meshgrid(np.linspace(0., 4*np.pi, 50),
                       np.linspace(.5 * np.pi, 4.5*np.pi, 50),
                       sparse=False, indexing='ij')
    error = np.zeros(x.shape)
    data = np.sin(x) + np.cos(y)
    data = IQazimuthal(intensity=data, error=error, qx=x, qy=y)

    plot_IQazimuthal(data, filename=filename, backend=backend)
    fileCheckAndRemove(filename, False)


def test_detector():
    pass


if __name__ == '__main__':
    pytest.main([__file__])
