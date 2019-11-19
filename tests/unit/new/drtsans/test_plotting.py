from drtsans.plots import plot_IQmod
from drtsans.dataobjects import IQmod
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


def test_IQazimuthal():
    pass


def test_detector():
    pass


if __name__ == '__main__':
    pytest.main([__file__])
