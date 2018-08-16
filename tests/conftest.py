from __future__ import (absolute_import, division, print_function)

import sys
import os
import pytest
from os.path import join as pjoin
from collections import namedtuple
import time
import mantid.simpleapi as mtds

# Resolve the path to the "external data"
this_module_path = sys.modules[__name__].__file__
parent_dir = pjoin(os.path.dirname(this_module_path), os.pardir)
data_dir = pjoin(parent_dir, 'data')


def timeit(a_function):
    """Calculate execution time"""
    def timed(*args, **kw):
        ts = time.time()
        result = a_function(*args, **kw)
        te = time.time()
        print('%r  %2.2f ms' % (a_function.__name__, (te - ts) * 1000))
        return result
    return timed


@pytest.fixture(scope='session')
def eqsans_f():
    return dict(data=pjoin(data_dir, 'eqsans', 'EQSANS_68168_event.nxs'),
                beamcenter=pjoin(data_dir, 'eqsans', 'EQSANS_68183_event.nxs'),
                darkcurrent=pjoin(data_dir, 'eqsans', 'EQSANS_68200_event.nxs')
                )


@pytest.fixture(scope='session')
def eqsans_w(eqsans_f):
    """Load EQSANS files into workspaces"""
    return {k: mtds.LoadEventNexus(v, OutputWorkspace=k)
            for (k, v) in eqsans_f.items()}


@pytest.fixture(scope='session')
def eqsans_p():
    """ Default parameters. Usually this comes from the parameters file """
    return dict(
        tubes_to_mask="1,48,53,54,85,123,130,137",
    )


@pytest.fixture(scope='session')
def porasil_slice1m():
    """EQSANS reduction benchmark. See porasil_slice1m.help

    File are loaded only once. For entries pointing to the same file, such as
    'dbc' and 'dbte', the file is loaded only once.

    Request to access a specific workspace through porasil_slice1m.w
    trigger cloning the workspace. Cloning is much quicker than reloading
    the file.

    """

    _help = """s: sample
    m: mask
    dc: dark current
    se: sensitivity
    dbc: direct_beam_center
    dbts: direct beam transmission sample
    dbte: direct beam transmission empty
    b: background
    bdbts: background direct beam transmission sample
    bdbte: background_direct_beam_transmission_empty
    """

    ipts = '/SNS/EQSANS/IPTS-20196'
    shared = '/SNS/EQSANS/shared/NeXusFiles/EQSANS'

    def fr(run_number):
        """Nexus file path from run number"""
        return pjoin(ipts, 'nexus', 'EQSANS_{}.nxs.h5'.format(run_number))

    # run numbers
    r = dict(s='92164',  # sample
             dc='89157',  # dark current
             dbc='92160',  # direct_beam_center
             dbts='92162',  # direct beam transmission sample
             dbte='92160',  # direct beam transmission empty
             b='92163',  # background
             bdbts='92161',  # background direct beam transmission sample
             bdbte='92160'  # background_direct_beam_transmission_empty
             )

    # Absolute path to benchmark files
    f = dict(s=fr('92164'),  # sample
             m=pjoin(shared, '2017B_mp/beamstop60_mask_4m.nxs'),  # mask
             dc=pjoin(shared, '2017B_mp/EQSANS_89157.nxs.h5'),  # dark current
             se=pjoin(shared, '2017B_mp/Sensitivity_patched_thinPMMA_1o3m_87680_event.nxs'),  # noqa: E501
             dbc=fr('92160'),  # direct_beam_center
             dbts=fr('92161'),  # direct beam transmission sample
             dbte=fr('92160'),  # direct beam transmission empty
             b=fr('92163'),  # background
             bdbts=fr('92161'),  # background direct beam transmission sample
             bdbte=fr('92160')  # background_direct_beam_transmission_empty
             )

    class GetWS(object):
        """Serves workspaces by cloning them. Prevents overwritting

        :param f: dictionary of filenames to load
        """

        def __init__(self, f):
            processed = dict()
            self._w = dict()
            for k, v in f.items():
                for other_k, other_v in processed.items():
                    if v == other_v:
                        self._w[k] = mtds.CloneWorkspace(self._w[other_k],
                                                         OutputWorkspace='_'+k)
                if k not in self._w:
                    self._w[k] = mtds.Load(v, OutputWorkspace='_'+k)
                    processed[k] = v

        def __len__(self):
            return len(self._w)

        def keys(self):
            return self._w.keys()

        def __getitem__(self, item):
            return mtds.CloneWorkspace(self._w[item],
                                       OutputWorkspace=item)

        def __setitem__(self, key, value):
            msg = "'GetWS' object does not support item assignment"
            raise TypeError(msg)

    ret_val = namedtuple('ret_val', 'ipts shared help r f w')
    return ret_val(ipts=ipts, shared=shared, r=r, f=f, w=GetWS(f), help=_help)
