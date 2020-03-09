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


def plot_IQmod(workspaces, filename, loglog=True, backend='d3',
               errorbar_kwargs={}, **kwargs):
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
    errorbar_kwargs: dict
        Optional arguments to :py:obj:`matplotlib.axes.Axes.errorbar`
        Can be a comma separated list for each workspace
        e.g. ``{'label':'main,wing,both', 'color':'r,b,g', 'marker':'o,v,.'}``
    kwargs: dict
        Additional key word arguments for :py:obj:`matplotlib.axes.Axes`

    '''
    backend = Backend.getMode(backend)
    for workspace in workspaces:
        datatype = getDataType(workspace)
        if datatype != DataType.IQ_MOD:
            raise RuntimeError('Do not know how to plot type="{}"'.format(datatype))

    fig, ax = plt.subplots()
    handles = []
    for n, workspace in enumerate(workspaces):
        eb, _, _ = ax.errorbar(workspace.mod_q, workspace.intensity, yerr=workspace.error)
        for key in errorbar_kwargs:
            value = [v.strip() for v in errorbar_kwargs[key].split(',')]
            plt.setp(eb, key, value[min(n, len(value)-1)])
    ax.set_xlabel(_q_label(backend))
    ax.set_ylabel('Intensity')
    if loglog:
        ax.set_xscale('log')
        ax.set_yscale('log')

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels)

    if kwargs:
        plt.setp(ax, **kwargs)

    _saveFile(fig, filename, backend)


def plot_IQazimuthal(workspace, filename, backend='d3',
                     imshow_kwargs={}, **kwargs):
    '''Save a plot representative of the supplied workspace

    Parameters
    ----------
    workspace: ~drtsans.dataobjects.IQazimuthal
        The workspace to plot. This assumes the data is binned on a constant grid.
    filename: str
        The name of the file to save to. For the :py:obj:`~Backend.MATPLOTLIB`
        backend, the type of file is determined from the file extension
    backend: Backend
        Which backend to save the file using
    imshow_kwargs: dict
        Optional arguments to :py:obj:`matplotlib.axes.Axes.imshow` e.g. ``{"norm": LogNorm()}``
    kwargs: dict
        Additional key word arguments for :py:obj:`matplotlib.axes.Axes`
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
                    origin='lower', aspect='auto', **imshow_kwargs)
    fig.colorbar(pcm, ax=ax)
    ax.set_xlabel(_q_label(backend, 'x'))
    ax.set_ylabel(_q_label(backend, 'y'))

    if kwargs:
        plt.setp(ax, **kwargs)

    _saveFile(fig, filename, backend)


def plot_detector(input_workspace, filename=None, backend='d3', axes_mode='tube-pixel',
                  panel_name=None, imshow_kwargs={'norm': LogNorm(vmin = 1)}):
    r"""
    Save a 2D plot representative of the supplied workspace

    Parameters
    ----------
    input_workspace: str, ~mantid.api.MatrixWorkspace
        The workspace to plot
    filename: str
        The name of the file to save to. For the :py:obj:`~Backend.MATPLOTLIB`
        backend, the type of file is determined from the file extension
    backend: Backend
        Which backend to save the file using
    axes_mode: str
        Plot intensities versus different axes. Options are: 'xy' for plotting versus pixel coordinates;
       'tube-pixel' for plotting versus tube and pixel index.
    panel_name: str
        Name of the double panel detector array. If :py:obj:`None`, plots will be generated for all arrays.
    imshow_kwargs: dict
        Optional arguments to matplotlib.axes.Axes.imshow
    """
    workspace = mtd[str(input_workspace)]
    backend = Backend.getMode(backend)
    detector_names = [panel_name, ] if panel_name is not None else panel_names(input_workspace)
    fig = plt.figure()
    for i_detector, detector_name in enumerate(detector_names):
        collection = TubeCollection(workspace, detector_name).sorted(view='decreasing X')
        data = np.sum(np.array([tube.readY for tube in collection]), axis=-1)  # sum intensities for each pixel
        if isinstance(imshow_kwargs.get('norm', None), LogNorm) is True:
            data[data < 1e-10] = 1e-10  # no negative values when doing a logarithm plot
        mask = np.array([tube.isMasked for tube in collection])
        data = np.ma.masked_where(mask, data)
        # Add subfigure
        axis = fig.add_subplot(len(detector_names), 1, i_detector + 1)
        if axes_mode == 'tube-pixel':
            image = axis.imshow(np.transpose(data), aspect='auto', origin='lower', **imshow_kwargs)
            axis_properties = {'set_xlabel': 'tube', 'set_ylabel': 'pixel', 'set_title': f'{detector_name}'}
        elif axes_mode == 'xy':
            n_pixels = len(collection[0])  # number of pixels in the first tube
            x = np.array([tube.position[0] * np.ones(n_pixels) for tube in collection])
            # BOTTLENECK-1
            y = np.array([tube.pixel_y for tube in collection])
            # BOTTLENECK-2 (but 6 times faster than BOTTLENECK-1)
            image = axis.pcolormesh(x, y, data)
            axis_properties = {'set_xlabel': 'X', 'set_ylabel': 'Y', 'set_title': f'{detector_name}'}
        else:
            raise ValueError('Unrecognized axes_mode. Valid options are "tube-pixel" and "xy"')
        image.cmap.set_bad(alpha=0.5)
        [getattr(axis, prop)(value) for prop, value in axis_properties.items()]
        fig.colorbar(image, ax=axis)
    fig.tight_layout()
    if filename is not None:
        _saveFile(fig, filename, backend)
