from __future__ import (absolute_import, division, print_function)

import sys
import os

import pytest
from os.path import join as pjoin
from collections import namedtuple
import mantid.simpleapi as mtds
from ornl.settings import amend_config

# Resolve the path to the "external data"
this_module_path = sys.modules[__name__].__file__
parent_dir = pjoin(os.path.dirname(this_module_path), os.pardir)
data_dir = pjoin(parent_dir, 'data')


def fr(ipts_number, run_number):
    """Nexus file path from run number"""
    return pjoin(ipts_number, 'nexus', 'EQSANS_{}.nxs.h5'.format(run_number))


ret_val = namedtuple('ret_val', 'ipts shared help r f w')


class GetWS(object):
    """Serves workspaces by cloning them. Prevents overwritting

    Parameters
    ----------
    f: dict
        Filenames to load
    p: str
        Prefix for workspace names. Prevents name overwriting of workspaces if
        this class is used in different fixtures
    loaders: dict
        Names of the Mantid algorithms used to load each file
    """

    def __init__(self, f, p, loaders=None):
        processed = dict()
        self._w = dict()
        for k, v in f.items():
            name = '_{}_{}'.format(p, k)  # workspace name. Begins with '_'
            for other_k, other_v in processed.items():
                if v == other_v:
                    self._w[k] = mtds.CloneWorkspace(self._w[other_k],
                                                     OutputWorkspace=name)
            if k not in self._w:
                loader_algm = mtds.Load
                if loaders is not None:
                    loader_algm = getattr(mtds, loaders[k])
                self._w[k] = loader_algm(v, OutputWorkspace=name)
                processed[k] = v

    def __len__(self):
        return len(self._w)

    def keys(self):
        return self._w.keys()

    def __getitem__(self, item):
        name = self._w[item].name()[1:]  # drop the intial '_'
        return mtds.CloneWorkspace(self._w[item],
                                   OutputWorkspace=name)

    def __setitem__(self, key, value):
        msg = "'GetWS' object does not support item assignment"
        raise TypeError(msg)


@pytest.fixture(scope='session')
def refd():
    """Directory locations for reference data

    Returns
    -------
    namedtuple
        refd.data, refd.legacy, refd.new,
        refd.legacy.biosans, refd.legacy.gpsans, refd.legacy.eqsans
        refd.new.biosans, refd.new.gpsans, refd.new.eqsans
        """
    d_leg = pjoin(data_dir, 'legacy', 'ornl', 'sans')
    d_new = pjoin(data_dir, 'new', 'ornl', 'sans')
    rett = namedtuple('rett', 'data legacy new')
    legt = namedtuple('legt', 'sans biosans gpsans eqsans')
    newt = namedtuple('newt', 'sans biosans gpsans eqsans')
    return rett(data_dir,
                legt(d_leg,
                     pjoin(d_leg, 'hfir', 'biosans'),
                     pjoin(d_leg, 'hfir', 'gpsans'),
                     pjoin(d_leg, 'sns', 'eqsans')),
                newt(d_new,
                     pjoin(d_new, 'hfir', 'biosans'),
                     pjoin(d_new, 'hfir', 'gpsans'),
                     pjoin(d_new, 'sns', 'eqsans')))


@pytest.fixture(scope='session')
def eqsans_f():
    return dict(data='EQSANS_68168',
                beamcenter='EQSANS_68183',
                darkcurrent='EQSANS_68200')


@pytest.fixture(scope='session')
def eqsans_w(eqsans_f):
    r"""Load EQSANS files into workspaces"""
    with amend_config({'instrumentName': 'EQSANS',
                       'datasearch.searcharchive': 'on'}):
        return {k: mtds.LoadEventNexus(v, OutputWorkspace=k)
                for (k, v) in eqsans_f.items()}


@pytest.fixture(scope='session')
def eqsans_p():
    """ Default parameters. Usually this comes from the parameters file """
    return dict(
        tubes_to_mask="1,48,53,54,85,123,130,137",
    )


@pytest.fixture(scope='session')
def biosans_f():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'biosans')
    return dict(
        beamcenter=pjoin(dd, 'BioSANS_exp402_scan0006_0001.xml'),
    )


@pytest.fixture(scope='session')
def gpsans_f():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'gpsans')
    return dict(
        beamcenter=pjoin(dd, 'CG2_exp325_scan0020_0001.xml'),
        beamcenter_off_setted=pjoin(dd, 'CG2_exp245_scan0007_0001.xml'),
        sample_transmission=pjoin(dd, 'CG2_exp245_scan0009_0001.xml'),
        sample_scattering=pjoin(dd, 'CG2_exp245_scan0010_0001.xml'),
        dark_current=pjoin(dd, 'CG2_exp244_scan0001_0001.xml'),
    )


@pytest.fixture(scope='session')
def gpsans_full_dataset():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'gpsans')
    return dict(
        sample_scattering_list=[
            pjoin(dd, 'CG2_exp245_scan0010_0001.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0002.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0003.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0004.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0005.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0006.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0007.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0008.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0009.xml'),
            pjoin(dd, 'CG2_exp245_scan0010_0010.xml'),
        ],
        background_scattering=pjoin(dd, 'CG2_exp245_scan0005_0001.xml'),
        sample_transmission=pjoin(dd, 'CG2_exp245_scan0009_0001.xml'),
        background_transmission=pjoin(dd, 'CG2_exp245_scan0004_0001.xml'),
        empty_transmission=pjoin(dd, 'CG2_exp245_scan0004_0001.xml'),
        beamcenter=pjoin(dd, 'CG2_exp245_scan0007_0001.xml'),
        dark_current=pjoin(dd, 'CG2_exp244_scan0001_0001.xml'),
    )


@pytest.fixture(scope='session')
<<<<<<< HEAD
def biosans_sensitivity_dataset():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'biosans')
    return dict(
        dark_current=pjoin(dd, 'BioSANS_exp327_scan0014_0001.xml'),
        flood=pjoin(dd, 'BioSANS_exp327_scan0066_0001.xml'),
        flood_beamcenter=pjoin(dd, 'BioSANS_exp327_scan0028_0001.xml'),
        empty_transmission=pjoin(dd, 'BioSANS_exp327_scan0028_0001.xml'),
        flood_mask=pjoin(dd, 'BioSANS_exp327_scan0066_0001_mask.xml'),
    )


@pytest.fixture(scope='session')
=======
>>>>>>> 1775da7... Added GPSANS sensitivity dataset
def gpsans_sensitivity_dataset():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'gpsans')
    return dict(
        dark_current=pjoin(dd, 'CG2_exp206_scan0038_0001.xml'),
<<<<<<< HEAD
        flood_trans=pjoin(dd, 'CG2_exp206_scan0017_0001.xml'),
        flood_trans_0_beamcenter=pjoin(dd, 'CG2_exp206_scan0016_0001.xml'),
        flood_trans_0_mask=pjoin(
            dd, 'CG2_exp206_scan0017_0001_mask_beamstop.xml'),
        flood_trans_200=pjoin(dd, 'CG2_exp206_scan0019_0001.xml'),
        flood_trans_200_beamcenter=pjoin(dd, 'CG2_exp206_scan0018_0001.xml'),
        flood_trans_200_mask=pjoin(
            dd, 'CG2_exp206_scan0019_0001_mask_beamstop.xml'),
        flood_trans_400=pjoin(dd, 'CG2_exp206_scan0021_0001.xml'),
        flood_trans_400_beamcenter=pjoin(dd, 'CG2_exp206_scan0020_0001.xml'),
        flood_trans_400_mask=pjoin(
            dd, 'CG2_exp206_scan0021_0001_mask_beamstop.xml'),
=======
        flood_trans_0=pjoin(dd, 'CG2_exp206_scan0017_0001.xml'),
        flood_trans_0_beamcenter=pjoin(dd, 'CG2_exp206_scan0016_0001.xml'),
        flood_trans_0_mask=pjoin(dd, 'CG2_exp206_scan0017_0001_mask.xml'),
        flood_trans_200=pjoin(dd, 'CG2_exp206_scan0019_0001.xml'),
        flood_trans_200_beamcenter=pjoin(dd, 'CG2_exp206_scan0018_0001.xml'),
        flood_trans_200_mask=pjoin(dd, 'CG2_exp206_scan0019_0001_mask.xml'),
        flood_trans_400=pjoin(dd, 'CG2_exp206_scan0021_0001.xml'),
        flood_trans_400_beamcenter=pjoin(dd, 'CG2_exp206_scan0020_0001.xml'),
        flood_trans_400_mask=pjoin(dd, 'CG2_exp206_scan0021_0.xml'),
>>>>>>> 1775da7... Added GPSANS sensitivity dataset
    )


@pytest.fixture(scope='session')
def frame_skipper():
    """Data and monitor with frame skipping
        """

    _help = """s: sample
    """

    ipts = '/SNS/EQSANS/IPTS-19150'
    shared = '/SNS/EQSANS/IPTS-19150/shared'

    # run numbers
    r = dict(s='92353',  # sample
             mo='92353',  # monitors
             )

    # Absolute path to benchmark files
    f = dict(s=fr(ipts, '92353'),  # sample
             mo=fr(ipts, '92353'),  # monitors
             )

    # Loader algorithms for the benchmark files
    lds = dict(s='Load',
               mo='LoadNexusMonitors',
               )

    return ret_val(ipts=ipts, shared=shared, r=r, f=f,
                   w=GetWS(f, 'frame_skipper', loaders=lds), help=_help)


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
    f = dict(s=fr(ipts, '92164'),  # sample
             m=pjoin(shared, '2017B_mp/beamstop60_mask_4m.nxs'),  # mask
             dc=pjoin(shared, '2017B_mp/EQSANS_89157.nxs.h5'),  # dark current
             se=pjoin(shared, '2017B_mp/Sensitivity_patched_thinPMMA_1o3m_87680_event.nxs'),  # noqa: E501
             dbc=fr(ipts, '92160'),  # direct_beam_center
             dbts=fr(ipts, '92161'),  # direct beam transmission sample
             dbte=fr(ipts, '92160'),  # direct beam transmission empty
             b=fr(ipts, '92163'),  # background
             bdbts=fr(ipts, '92161'),  # noqa: E501 background direct beam transmission sample
             bdbte=fr(ipts, '92160')  # noqa: E501 background_direct_beam_transmission_empty
             )

    lds = dict(s='Load',  # sample
               m='Load',  # mask
               dc='Load',  # dark current
               se='Load',  # sensitivity
               dbc='Load',  # direct_beam_center
               dbts='Load',  # direct beam transmission sample
               dbte='Load',  # direct beam transmission empty
               b='Load',  # background
               bdbts='Load',  # background direct beam transmission sample
               bdbte='Load'  # background_direct_beam_transmission_empty
               )

    return ret_val(ipts=ipts, shared=shared, r=r, f=f,
                   w=GetWS(f, 'porasil_slice1m', loaders=lds), help=_help)
