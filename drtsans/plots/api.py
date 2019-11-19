from drtsans.dataobjects import DataType, getDataType
from enum import Enum
import json
import numpy as np
import mpld3 # noqa E402
from mpld3 import plugins # noqa E402
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt # noqa E402


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
    MPLD3 = 'd3'     # read-only
    MATPLOTLIB = 'mpl'  # read and write

    def __str__(self):
        return self.value

    @staticmethod
    def getMode(mode):
        '''Private function to convert anything into :py:obj:`HidraProjectFileMode`'''
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
    if backend == Backend.MATPLOTLIB:
        figure.savefig(filename)
        if show:
            raise RuntimeError('Cannot show data with matplotlib backend')
    else:
        if not filename.endswith('json'):
            raise RuntimeError('File "{}" must have ".json" suffix'.format(filename))

        plugins.connect(figure, plugins.MousePosition(fontsize=14, fmt=".0f"))
        with open(filename, 'w') as outfile:
            mpld3.save_json(figure, outfile)

        if show:
            mpld3.show(figure)


def plot_IQmod(workspaces, filename='', backend='d3'):
    backend = Backend.getMode(backend)
    for workspace in workspaces:
        datatype = getDataType(workspace)
        if datatype != DataType.IQ_MOD:
            raise RuntimeError('Do not know how to plot type="{}"'.format(datatype))

    fig, ax = plt.subplots()
    for workspace in workspaces:
        ax.errorbar(workspace.mod_q, workspace.intensity, yerr=workspace.error)
    ax.set_xlabel(r'$Q (\AA^{-1})$')
    ax.set_ylabel('Intensity')

    _saveFile(fig, filename, backend)


def plot_IQazimuthal():
    pass


def plot_detector():
    pass
