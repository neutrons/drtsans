from drtsans.plots import plot_IQmod, plot_IQazimuthal, plot_detector
from drtsans.dataobjects import IQmod, IQazimuthal
from mantid.simpleapi import LoadEmptyInstrument, LoadNexus
import numpy as np
import os
import pytest
import matplotlib.pyplot as plt


def fileCheckAndRemove(filename, remove=True):
    '''Convienience function for doing simple checs that the file was created.
    The ``remove`` option is available to make debugging new tests easier.'''
    assert os.path.exists(filename), '"{}" does not exist'.format(filename)
    assert os.path.getsize(filename) > 100, '"{}" is too small'.format(filename)
    if remove:
        os.remove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQmod.png'),
                          ('d3', 'test_IQmod.json')],
                         ids=['mpl', 'd3'])
def test_IQmod(backend, filename):
    '''Test plotting single a IQmod dataset'''
    x = np.linspace(0., 4*np.pi, 50)
    e = np.zeros(50) + .1
    data = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    plot_IQmod([data], filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQmod_multi.png'),
                          ('d3', 'test_IQmod_multi.json')],
                         ids=['mpl', 'd3'])
def test_IQmod_multi(backend, filename):
    '''Test over-plotting multiple IQmod datasets'''
    x = np.linspace(0., 4*np.pi, 50)
    e = np.zeros(50) + .1
    data1 = IQmod(intensity=np.sin(x), error=e, mod_q=x)
    data2 = IQmod(intensity=np.cos(x), error=e, mod_q=x)

    plot_IQmod([data1, data2], filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQazimuthal_1d.png'),
                          ('d3', 'test_IQazimuthal_1d.json')],
                         ids=['mpl', 'd3'])
def test_IQazimuthal_1d(backend, filename):
    '''Test plotting IQazimuthal with 1d Qx and Qy arrays'''
    x = np.linspace(0., 4*np.pi, 50)
    y = np.linspace(.5 * np.pi, 4.5*np.pi, 50)
    error = np.zeros((50, 50))
    data = error[:]
    for i in range(x.size):
        for j in range(y.size):
            data[i, j] = np.sin(x[i]) + np.cos(y[j])
    data = IQazimuthal(intensity=data, error=error, qx=x, qy=y)

    plot_IQazimuthal(data, filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQazimuthal_2d.png'),
                          ('d3', 'test_IQazimuthal_2d.json')],
                         ids=['mpl', 'd3'])
def test_IQazimuthal_2d(backend, filename):
    '''Test plotting IQazimuthal with 2d Qx and Qy arrays'''
    x, y = np.meshgrid(np.linspace(0., 4*np.pi, 50),
                       np.linspace(.5 * np.pi, 4.5*np.pi, 50),
                       sparse=False, indexing='ij')
    error = np.zeros(x.shape)
    data = np.sin(x) + np.cos(y)
    data = IQazimuthal(intensity=data, error=error, qx=x, qy=y)

    plot_IQazimuthal(data, filename=filename, backend=backend)
    plt.close()
    fileCheckAndRemove(filename, False)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_IQazimuthal_2d_selections.png'),
                          ('d3', 'test_IQazimuthal_2d_selections.json')],
                         ids=['mpl', 'd3'])
def test_IQazimuthal_2d_selections(backend, filename):
    '''Test plotting IQazimuthal with 2d Qx and Qy arrays'''
    x, y = np.meshgrid(np.linspace(0., 4*np.pi, 50),
                       np.linspace(.5 * np.pi, 4.5*np.pi, 50),
                       sparse=False, indexing='ij')
    error = np.zeros(x.shape)
    data = np.sin(x) + np.cos(y)
    data = IQazimuthal(intensity=data, error=error, qx=x, qy=y)

    plot_IQazimuthal(data, filename=filename, backend=backend,
                     qmin=0., qmax=2., wedges=((-45., 45.),))
    plt.close()
    fileCheckAndRemove(filename, False)


@pytest.mark.parametrize('backend, filename',
                         [('mpl', 'test_detector.png'),
                          ('d3', 'test_detector.json')],
                         ids=['mpl', 'd3'])
def test_detector(backend, filename):
    '''Test plotting in detector space from a mantid workspace'''
    workspace = LoadEmptyInstrument(InstrumentName='CG3')  # this will load monitors as well
    plot_detector(workspace, filename, backend)
    plt.close()
    fileCheckAndRemove(filename, False)


def test_xaxis_direction(reference_dir):
    r"""Test values of X-axis in plot_detector decrease when looking at the picture from left to right"""
    # wing_detector.nxs contains intensities for the wing detector that can be plotted
    workspace = LoadNexus(os.path.join(reference_dir.new.sans, 'plots', 'wing_detector.nxs'))
    filename = 'test_xaxis_direction.png'
    plot_detector(workspace, filename=filename, backend='mpl', panel_name='wing_detector', axes_mode='xy')
    plt.close()
    fileCheckAndRemove(filename, remove=True)


if __name__ == '__main__':
    pytest.main([__file__])
