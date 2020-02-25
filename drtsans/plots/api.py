from enum import Enum
import json
import numpy as np
import matplotlib
import warnings
warnings.simplefilter('ignore', UserWarning)
matplotlib.use('Agg')
import matplotlib.pyplot as plt # noqa E402
from matplotlib.colors import LogNorm # noqa E402
import mpld3 # noqa E402
from mpld3 import plugins # noqa E402

from mantid.api import mtd # noqa E402
from drtsans.tubecollection import TubeCollection # noqa E402
from drtsans.dataobjects import DataType, getDataType # noqa E402
from drtsans.geometry import panel_names # noqa E402


__all__ = ['plot_IQmod', 'plot_IQazimuthal', 'plot_detector']


# mpld3 taken from hack from https://github.com/mpld3/mpld3/issues/434#issuecomment-381964119
if mpld3.__version__ == '0.3':
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)
    from mpld3 import _display  # noqa E402
    _display.NumpyEncoder = NumpyEncoder


class Backend(Enum):
    '''Class for denoting which back-end to save the plot using'''
    MPLD3 = 'd3'     # read-only
    MATPLOTLIB = 'mpl'  # read and write

    def __str__(self):
        return self.value

    @staticmethod
    def getMode(mode):
        '''Public function to convert anything into :py:obj:`Backend`'''
        if mode in Backend:
            return mode
        else:
            mode = str(mode)
            if mode == 'mpl' or mode == 'matplotlib':
                return Backend.MATPLOTLIB
            elif mode == 'mpld3' or mode == 'd3':
                return Backend.MPLD3
            else:
                return Backend[mode.upper()]


def _saveFile(figure, filename, backend, show=False):
    '''Convenience method for the common bits of saving the file based on
    the selected backend.

    Parameters
    ----------
    figure: ~matplotlib.pyplot.figure
        The figure to save to a file
    filename: str
        The name of the file to save to
    backend: Backend
        Which :py:obj:`Backend` to use for saving
    show: bool
        Whether or not to show the figure rather than saving. This is only
        available if the :py:obj:`~Backend.MPLD3` backend is selected.
    '''
    if backend == Backend.MATPLOTLIB:
        if filename:
            figure.savefig(filename)
        if show:
            if 'inline' in matplotlib.get_backend():  # ipython notebook
                figure.show()
            else:
                raise RuntimeError('Cannot show data with matplotlib backend')
    else:
        if not filename.endswith('json'):
            raise RuntimeError('File "{}" must have ".json" suffix'.format(filename))

        plugins.connect(figure, plugins.MousePosition(fontsize=14, fmt=".0f"))
        with open(filename, 'w') as outfile:
            mpld3.save_json(figure, outfile)

        if show:
            mpld3.show(figure)


def _q_label(backend, subscript=''):
    '''mpld3 doesn't currently support latex markup. This generates an
    acceptable label for the supplied backend.

    Parameters
    ----------
    backend: Backend
        Which :py:obj:`Backend` to generate the caption for
    subscript: str
        The subscript on the "Q" label. If none is specified then no
        subscript will be displayed
    '''
    label = 'Q'
    if subscript:
        label += '_' + str(subscript)

    if backend == Backend.MATPLOTLIB:
        return '$' + label + r' (\AA^{-1})$'
    else:  # mpld3
        return label + ' (1/{})'.format(u'\u212B')


def plot_IQmod(workspaces, filename, loglog=True, backend='d3'):
    '''Save a plot representative of the supplied workspaces

    Parameters
    ----------
    workspaces: list, tuple
        A collection of :py:obj:`~drtsans.dataobjects.IQmod` workspaces to
        plot. If only one is desired, it must still be supplied in a
        :py:obj:`list` or :py:obj:`tuple`.
    filename: str
        The name of the file to save to. For the :py:obj:`~Backend.MATPLOTLIB`
        backend, the type of file is determined from the file extension
    loglog: bool
        If true will set both axis to logarithmic, otherwise leave them as linear
    backend: Backend
        Which backend to save the file using
    '''
    backend = Backend.getMode(backend)
    for workspace in workspaces:
        datatype = getDataType(workspace)
        if datatype != DataType.IQ_MOD:
            raise RuntimeError('Do not know how to plot type="{}"'.format(datatype))

    fig, ax = plt.subplots()
    for workspace in workspaces:
        ax.errorbar(workspace.mod_q, workspace.intensity, yerr=workspace.error)
    ax.set_xlabel(_q_label(backend))
    ax.set_ylabel('Intensity')
    if loglog:
        ax.set_xscale('log')
        ax.set_yscale('log')

    _saveFile(fig, filename, backend)


def plot_IQazimuthal(workspace, filename, backend='d3'):
    '''Save a plot representative of the supplied workspace

    Parameters
    ----------
    workspaces: ~drtsans.dataobjects.IQazimuthal
        The workspace to plot. This assumes the data is binned on a constant grid.
    filename: str
        The name of the file to save to. For the :py:obj:`~Backend.MATPLOTLIB`
        backend, the type of file is determined from the file extension
    backend: Backend
        Which backend to save the file using
    '''
    backend = Backend.getMode(backend)
    datatype = getDataType(workspace)
    if datatype != DataType.IQ_AZIMUTHAL:
        raise RuntimeError('Do not know how to plot type="{}"'.format(datatype))

    qxmin = workspace.qx.min()
    qxmax = workspace.qx.max()
    qymin = workspace.qy.min()
    qymax = workspace.qy.max()

    fig, ax = plt.subplots()
    current_cmap = matplotlib.cm.get_cmap()
    current_cmap.set_bad(color='grey')
    pcm = ax.imshow(workspace.intensity.T, extent=(qxmin, qxmax, qymin, qymax),
                    norm=LogNorm(), origin='lower', aspect='auto')
    fig.colorbar(pcm, ax=ax)
    ax.set_xlabel(_q_label(backend, 'x'))
    ax.set_ylabel(_q_label(backend, 'y'))

    _saveFile(fig, filename, backend)


def plot_detector(input_workspace, filename=None, backend='d3',
                  imshow_kwargs={'norm': LogNorm(vmin=1)}):
    r"""Save a plot representative of the supplied workspace

     Parameters
     ----------
     workspaces: str, ~mantid.api.MatrixWorkspace
         The workspace to plot
     filename: str
         The name of the file to save to. For the :py:obj:`~Backend.MATPLOTLIB`
         backend, the type of file is determined from the file extension
     backend: Backend
         Which backend to save the file using
     imshow_kwargs: dict
         Optional arguments to matplotlib.axes.Axes.imshow
     """
    workspace = mtd[str(input_workspace)]
    backend = Backend.getMode(backend)
    detector_names = panel_names(input_workspace)
    fig = plt.figure()
    for i_detector, detector_name in enumerate(detector_names):
        collection = TubeCollection(workspace, detector_name).sorted(view='decreasing X')
        data = np.sum(np.array([tube.readY for tube in collection]), axis=-1)  # sum intensities for each pixel
        data[data < 1e-10] = 1e-10
        mask = np.array([tube.isMasked for tube in collection])
        data = np.transpose(np.ma.masked_where(mask, data))
        # Add subfigure
        axis = fig.add_subplot(len(detector_names), 1, i_detector + 1)
        image = axis.imshow(data, aspect='auto', origin='lower', **imshow_kwargs)
        image.cmap.set_bad(alpha=0.5)
        axis_properties = {'set_xlabel': 'tube', 'set_ylabel': 'pixel', 'set_title': f'{detector_name}'}
        [getattr(axis, prop)(value) for prop, value in axis_properties.items()]
        fig.colorbar(image, ax=axis)
    fig.tight_layout()
    if filename is not None:
        _saveFile(fig, filename, backend)
