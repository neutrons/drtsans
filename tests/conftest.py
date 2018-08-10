from __future__ import (absolute_import, division, print_function)

import sys
import os
import pytest
from os.path import join as pjoin
import mantid.simpleapi as mtds
import time

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
    """EQSANS reduction benchmark"""

    def _porasil_slice1m(fs=None):
        """
        :param fs: load files into workspaces. None will not load anything,
        'all' will load all files, 'dc' will load the dark current,
        ('s', 'dc') will load the sample and the dark_current files
        """
        ipts = '/SNS/EQSANS/IPTS-20196'
        shared = '/SNS/EQSANS/shared/NeXusFiles/EQSANS'

        def fr(run_number):
            """Nexus file path from run number"""
            return pjoin(ipts, 'nexus', 'EQSANS_92164.nxs.h5')

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
                 se=pjoin(shared, '2017B_mp/Sensitivity_patched_thinPMMA_1o3m_87680_event.nxs'),  # sensitivity
                 dbc=fr('92160'),  # direct_beam_center
                 dbts=fr('92161'),  # direct beam transmission sample
                 dbte=fr('92160'),  # direct beam transmission empty
                 b=fr('92163'),  # background
                 bdbts=fr('92161'),  # background direct beam transmission sample
                 bdbte=fr('92160')  # background_direct_beam_transmission_empty
                 )
        # Workspaces
        w = dict()
        if fs is not None:
            if fs == 'all':
                w = {k: mtds.Load(v, OutputWorkspace=k) for (k, v) in f.items()}
            elif fs in f.keys():
                w[fs] = mtds.Load(f[fs], OutputWorkspace=fs)
            else:
                # fs is an iterable containing keys of dictionary 'f'
                w = {k: mtds.Load(f[k], OutputWorkspace=k) for k in fs}

        return dict(ipts=ipts, shared=shared, r=r, f=f, w=w)

    return _porasil_slice1m