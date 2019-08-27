from __future__ import (absolute_import, division, print_function)

import sys
import os

import pytest
import random
import string
import numpy as np
from os.path import join as pjoin
from collections import namedtuple
import mantid.simpleapi as mtds
from mantid.simpleapi import CreateWorkspace, LoadInstrument
from ornl.settings import amend_config

# Resolve the path to the "external data"
this_module_path = sys.modules[__name__].__file__
parent_dir = pjoin(os.path.dirname(this_module_path), os.pardir)

data_dir = '/SNS/EQSANS/shared/sans-backend/data'


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


@pytest.fixture(scope='module')
def cleanfile():
    '''Fixture that deletes registered files when the .py file is finished. It
    will cleanup on exception and will safely skip over files that do not
    exist.

    Usage:

    def test_something(cleanfile):
        cleanfile('/some/file/the/test.creates')
        # do stuff
    '''
    filenames = []

    def _cleanfile(filename):
        filenames.append(filename)
        return filename

    yield _cleanfile

    for name in filenames:
        if os.path.exists(name):
            os.unlink(name)


@pytest.fixture(scope='session')
def reference_dir():
    """A namedtuple with the directory **absolute** paths for test data

    Examples:
        reference_dir.data, topmost data directory data/
        reference_dir.legacy, data/legacy/ornl/sans/
        reference_dir.new, data/new/ornl/sans/
        reference_dir.legacy.biosans, reference_dir.legacy.gpsans, reference_dir.legacy.eqsans, are
            data/legacy/ornl/sans/hfir/biosans and so on.
        reference_dir.new.biosans, reference_dir.new.gpsans, reference_dir.new.eqsans, are
            data/new/ornl/sans/hfir/biosans and so on.

    Returns
    -------
    namedtuple
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
def eqsans_f(reference_dir):
    return dict(data=pjoin(reference_dir.new.eqsans, 'EQSANS_68168_event.nxs'),
                beamcenter=pjoin(reference_dir.new.eqsans, 'EQSANS_68183_event.nxs'),
                darkcurrent=pjoin(reference_dir.new.eqsans, 'EQSANS_68200_event.nxs'))


@pytest.fixture(scope='session')
def eqsans_w(reference_dir, eqsans_f):
    r"""Load EQSANS files into workspaces"""
    with amend_config(data_dir=reference_dir.new.eqsans):
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
        anisotropic=pjoin(dd, 'BioSANS_exp440_scan0022_0006.xml'),
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
        anisotropic=pjoin(dd, 'CG2_exp296_scan0166_0001.xml'),
        sample_scattering_2=pjoin(dd, 'CG2_exp325_scan0007_0001.xml'),
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
def gpsans_sensitivity_dataset():
    dd = pjoin(data_dir, 'new', 'ornl', 'sans', 'hfir', 'gpsans')
    return dict(
        dark_current=pjoin(dd, 'CG2_exp206_scan0038_0001.xml'),
        flood_trans_0=pjoin(dd, 'CG2_exp206_scan0017_0001.xml'),
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
    )


@pytest.fixture(scope='session')
def frame_skipperF(reference_dir):
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
    f = dict(s=pjoin(reference_dir.new.eqsans, 'EQSANS_92353.nxs.h5'),  # sample
             mo=pjoin(reference_dir.new.eqsans, 'EQSANS_92353.nxs.h5')  # monitors
             )

    # Loader algorithms for the benchmark files
    lds = dict(s='Load',
               mo='LoadNexusMonitors',
               )

    return ret_val(ipts=ipts, shared=shared, r=r, f=f,
                   w=GetWS(f, 'frame_skipper', loaders=lds), help=_help)


@pytest.fixture(scope='session')
def porasil_slice1m(reference_dir):
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
    f = dict(s=pjoin(reference_dir.new.eqsans, 'EQSANS_92164.nxs.h5'),  # sample
             m=pjoin(reference_dir.new.eqsans, '2017B_mp/beamstop60_mask_4m.nxs'),  # noqa: E501 mask
             dc=pjoin(reference_dir.new.eqsans, 'EQSANS_89157.nxs.h5'),  # dark current
             se=pjoin(reference_dir.new.eqsans, 'Sensitivity_patched_thinPMMA_1o3m_87680_event.nxs'),  # noqa: E501
             dbc=pjoin(reference_dir.new.eqsans, 'EQSANS_92160.nxs.h5'),  # noqa: E501 direct_beam_center
             dbts=pjoin(reference_dir.new.eqsans, 'EQSANS_92161.nxs.h5'),  # noqa: E501 direct beam transmission sample
             dbte=pjoin(reference_dir.new.eqsans, 'EQSANS_92160.nxs.h5'),  # noqa: E501 direct beam transmission empty
             b=pjoin(reference_dir.new.eqsans, 'EQSANS_92163.nxs.h5'),  # background
             bdbts=pjoin(reference_dir.new.eqsans, 'EQSANS_92161.nxs.h5'),  # noqa: E501 background direct beam transmission sample
             bdbte=pjoin(reference_dir.new.eqsans, 'EQSANS_92160.nxs.h5')  # noqa: E501 background_direct_beam_transmission_empty
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


@pytest.fixture()
def generic_IDF(request):
    '''
    generate a test IDF with a rectangular detector
    with Nx X Ny pixels

    Parameters
    ----------

    request is a dictionary containing the following keys:

        name: Name of the instrument     (default: GenericSANS)
        Nx : number of columns                      (default 3)
        Ny : number of rows                         (default 3)
        dx : width of a column in meters            (default 1)
        dy : height of a row in meters              (default 1)
        xc : distance of center along the x axis    (default 0)
        yc : distance of center along the y axis    (default 0)
        zc : distance of center along the z axis    (default 5)
        l1 : distance from source to sample       (default -11)

    Note that we use Mantid convention for the orientation
    '''
    # try to get the parent in case of sub-requests
    try:
        req_params = request.param
    except AttributeError:
        try:
            req_params = request._parent_request.param
        except AttributeError:
            req_params = dict()

    # use hidden attibutes to get data dimension, Nx and Ny can override this
    intensity = req_params.get('intensities', None)
    if intensity:
        try:
            Nx, Ny = intensity.shape[:2]  # numpy array
        except AttributeError:
            Nx = len(intensity)
            Ny = len(intensity[0])
    else:
        Nx, Ny = 3, 3

    # get the parameters from the request object
    params = {'name': req_params.get('name', 'GenericSANS'),
              'l1': -1. * abs(float(req_params.get('l1', -11.))),
              'Nx': int(req_params.get('Nx', Nx)),
              'Ny': int(req_params.get('Ny', Ny)),
              'dx': float(req_params.get('dx', 1.)),
              'dy': float(req_params.get('dy', 1.)),
              'xcenter': float(req_params.get('xc', 0.)),
              'ycenter': float(req_params.get('yc', 0.)),
              'zcenter': float(req_params.get('zc', 5.))}
    params['dx_mm'] = params['dx'] * 1000.
    params['dy_mm'] = params['dy'] * 1000.

    # check that nothing is crazy
    assert (params['Nx'] > 1 and params['Nx'] < 300)
    assert (params['Ny'] > 1 and params['Ny'] < 300)
    assert params['dx'] > 0.
    assert params['dy'] > 0.
    assert params['zcenter'] > 0.

    # derived parameters
    params['half_dx'] = params['dx'] * .5
    params['half_dy'] = params['dy'] * .5
    params['xstart'] = -(params['Nx']-1) * params['half_dx']
    params['ystart'] = -(params['Ny']-1) * params['half_dy']

    template_xml = '''<?xml version='1.0' encoding='UTF-8'?>
<instrument name="{name}" valid-from   ="1900-01-31 23:59:59"
                               valid-to     ="2100-12-31 23:59:59"
                               last-modified="2019-07-12 00:00:00">
    <!--DEFAULTS-->
    <defaults>
        <length unit="metre"/>
        <angle unit="degree"/>
        <reference-frame>
        <along-beam axis="z"/>
        <pointing-up axis="y"/>
        <handedness val="right"/>
        <theta-sign axis="x"/>
        </reference-frame>
    </defaults>

    <!--SOURCE-->
    <component type="moderator">
        <location z="{l1}"/>
    </component>
    <type name="moderator" is="Source"/>

    <!--SAMPLE-->
    <component type="sample-position">
        <location y="0.0" x="0.0" z="0.0"/>
    </component>
    <type name="sample-position" is="SamplePos"/>

    <!--RectangularDetector-->
    <component type="panel" idstart="0" idfillbyfirst="y" idstepbyrow="{Ny}">
        <location x="{xcenter}" y="{ycenter}" z="{zcenter}"
            name="detector1"
            rot="0.0" axis-x="0" axis-y="1" axis-z="0">
        </location>
    </component>

    <!-- Rectangular Detector Panel -->
    <type name="panel" is="rectangular_detector" type="pixel"
        xpixels="{Nx}" xstart="{xstart}" xstep="+{dx}"
        ypixels="{Ny}" ystart="{ystart}" ystep="+{dy}" >
        <properties/>
    </type>

    <!-- Pixel for Detectors-->
    <type is="detector" name="pixel">
        <cuboid id="pixel-shape">
            <left-front-bottom-point y="-{half_dy}" x="-{half_dx}" z="0.0"/>
            <left-front-top-point y="{half_dy}" x="-{half_dx}" z="0.0"/>
            <left-back-bottom-point y="-{half_dy}" x="-{half_dx}" z="-0.0001"/>
            <right-front-bottom-point y="-{half_dy}" x="{half_dx}" z="0.0"/>
        </cuboid>
        <algebra val="pixel-shape"/>
    </type>

    <parameter name="x-pixel-size">
        <value val="{dx_mm}"/>
    </parameter>

    <parameter name="y-pixel-size">
        <value val="{dy_mm}"/>
    </parameter>
</instrument>'''

    # return the completed template
    return template_xml.format(**params)


@pytest.fixture()
def generic_workspace(generic_IDF, request):
    '''
    generate a test IDF with a rectangular detector
    with Nx X Ny pixels

    Parameters
    ----------

    request is a dictionary containing the following keys:

        name: Name of the workspace and instrument
                                         (default: GenericSANS)
        axis_values : ndarray or 2d-list of the independent axis for the
             data. It will be copied across all spectra if only specified
             for one.               (default 0 for all spectra)
        intensities : ndarray or 2d/3d list of intensities for the
             instrument. Detector dimensions are inferred from the
             dimensionality. This will be linearized using `numpy.ravel`.
                          (default: zeros of dimension Nx x Ny)
        uncertainties : ndarray or 2d/3d list of intensities for the
             instrument. This will be linearized using `numpy.ravel`.
                                  (default: sqrt(intensities),
                                   or one if intensity is zero)
        axis_units : units for the independent axis
                                           (default wavelength)
        Nx : number of columns                      (default 3)
        Ny : number of rows                         (default 3)
        dx : width of a column in meters            (default 1)
        dy : height of a row in meters              (default 1)
        xc : distance of center along the x axis    (default 0)
        yc : distance of center along the y axis    (default 0)
        zc : distance of center along the z axis    (default 5)
        l1 : distance from source to sample       (default -11)

    Example
    -------
    For a workspace specified with the parameters (all other parameters are default)

    .. code-block:: python
       {'axis_values':[42.],
        'intensities': [[1.,4.],[9.,16.],[25.,36.]]}

    The intensities will be in 2x3 grid

    .. code-block:: python
       print(wksp.extractY().reshape(3,2)

    .. code-block::
       [[ 1.  4.]
       [ 9. 16.]
       [25. 36.]]

    which is vertically upside-down from how the data is on the detectors.
    The positions (in x,y) of the pixels (parallel to the previous array)

    .. code-block::
                y=-1.   y=0.   y=1.
       x=-0.5  id=0    id=2   id=4
       x= 0.5  id=1    id=3   id=5

    All z-values are 5.
    In the case of time-of-flight data, add more to `axis_values` and `axis_units='tof'`.

    Note that we use Mantid convention for the orientation
    '''
    try:
        req_params = request.param
    except AttributeError:
        try:
            req_params = request._parent_request.param
        except AttributeError:
            req_params = dict()

    name = req_params.get('name', 'GenericSANS')  # output workspace
    units = req_params.get('axis_units', 'wavelength')

    # get the supplied data
    x = req_params.get('axis_values', None)
    y = req_params.get('intensities', None)
    e = req_params.get('uncertainties', None)
    if y:
        try:
            Nx, Ny = y.shape[:2]
        except AttributeError:
            Nx = len(y)
            Ny = len(y[0])
            y = np.array(y)
    else:
        Nx = req_params.get('Nx', 3)  # must match generic_IDF
        Ny = req_params.get('Ny', 3)  # must match generic_IDF
        y = np.zeros((Nx, Ny), dtype=float)
    y = y.ravel()
    if e:
        e = np.array(e).ravel()
    else:
        e = np.sqrt(y)
        e[e == 0.] = 1.  # the default SANS likes
    if x:
        x = np.array(x).ravel()
    else:
        x = np.zeros(Nx * Ny, dtype=float)

    wksp = CreateWorkspace(DataX=x,
                           DataY=y,
                           DataE=e,
                           Nspec=Nx * Ny,
                           UnitX=units,
                           OutputWorkspace=name)
    LoadInstrument(Workspace=wksp, InstrumentXML=generic_IDF,
                   RewriteSpectraMap=True, InstrumentName=name)

    return wksp


@pytest.fixture(scope='session')
def serve_events_workspace(reference_dir):
    r"""
    Load an events workspace and cache it for future requests.

    If the same run is requested, the fixture clones the cached workspace,
    thus avoiding reloading the file.

    Parameters
    ----------
    run: str
        Instrument plus run number string, e.g 'EQSANS_92353'
    dd: str
        directory location where to find the file. Unnecessary if /SNS mounted.

    Returns
    -------
    EventsWorkspace
    """
    def wrapper(run, dd=reference_dir.new.eqsans):
        cache = wrapper._cache
        names = wrapper._names

        def uwd():
            while True:
                name = '__' + ''.join(random.choice(string.ascii_lowercase)
                                      for _ in range(9))
                if name not in names:
                    return name

        if cache.get(run, None) is None:
            with amend_config(data_dir=dd):
                cache[run] = mtds.LoadEventNexus(run, OutputWorkspace=uwd())
                names.append(cache[run])
        clone = mtds.CloneWorkspace(cache[run], OutputWorkspace=uwd())
        names.append(clone.name())
        return clone

    wrapper._cache = dict()  # caches the loaded original ws
    wrapper._names = list()  # stores names for all ws produced
    yield wrapper
    [mtds.DeleteWorkspace(name) for name in wrapper._names]
