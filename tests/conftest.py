import sys
import os
import re
import pytest
import random
import string
import numpy as np
from os.path import join as pjoin
from collections import namedtuple
import mantid.simpleapi as mtds
from mantid.simpleapi import CompareWorkspaces, CreateWorkspace, LoadInstrument, DeleteWorkspace
from drtsans.dataobjects import DataType, getDataType
from drtsans.settings import amend_config, unique_workspace_dundername

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
def gpsans_sensitivity_dataset(reference_dir):
    dd = reference_dir.new.gpsans
    return dict(
        dark_current=pjoin(dd, 'CG2_exp206_scan0038_0001.xml'),
        flood_trans_0=pjoin(dd, 'CG2_exp206_scan0017_0001.xml'),
        flood_trans_0_beamcenter=pjoin(dd, 'CG2_exp206_scan0016_0001.xml'),
        flood_trans_0_mask=pjoin(dd, 'CG2_exp206_scan0017_0001_mask_beamstop.xml'),
        flood_trans_200=pjoin(dd, 'CG2_exp206_scan0019_0001.xml'),
        flood_trans_200_beamcenter=pjoin(dd, 'CG2_exp206_scan0018_0001.xml'),
        flood_trans_200_mask=pjoin(dd, 'CG2_exp206_scan0019_0001_mask_beamstop.xml'),
        flood_trans_400=pjoin(dd, 'CG2_exp206_scan0021_0001.xml'),
        flood_trans_400_beamcenter=pjoin(dd, 'CG2_exp206_scan0020_0001.xml'),
        flood_trans_400_mask=pjoin(dd, 'CG2_exp206_scan0021_0001_mask_beamstop.xml'),
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


def _getDataDimensions(req_params):
    '''
    Determine the dimensionality of the data for use in generic_IDF and generic_workspace.
    The basic rules are that the dimensionality is taken from ``Nx`` and ``Ny``, then
    determined from ``intensities``, then is ``(3,3)``. If both dimensions and intensities
    are specified, dimensions wins, but it must be compatible with the number of values in
    intensities.

    Parameters
    ----------

    request is a dictionary containing the following keys:

        intensities : ndarray or 2d/3d list of intensities for the
             instrument. Detector dimensions are inferred from the
             dimensionality.
        Nx : number of columns                      (default 3)
        Ny : number of rows                         (default 3)

    Returns
    -------
    tuple
    Nx, Ny
    '''
    intensity = req_params.get('intensities', None)
    if intensity is None:
        Nx = int(req_params.get('Nx', 3))
        Ny = int(req_params.get('Ny', 3))
        return Nx, Ny, 1
    else:
        # force it to be a numpy array
        # this is a no-op if it is already the right type
        intensity = np.array(intensity)

        Nx = req_params.get('Nx', None)
        Ny = req_params.get('Ny', None)
        if (Nx is not None) and (Ny is not None):
            Nx = int(Nx)
            Ny = int(Ny)
            if intensity.size % (Nx * Ny) != 0:
                raise RuntimeError('Supplied Nx={}, Ny={} not compatible with '
                                   'intensities[{}]'.format(Nx, Ny, intensity.shape))
            else:
                return Nx, Ny, int(intensity.size / (Nx * Ny))
        else:
            if len(intensity.shape) == 3:
                return intensity.shape
            else:
                Nx, Ny = intensity.shape[:2]  # Nx, Ny
                return Nx, Ny, 1


def idf_xml_factory(idf_xml_name, request):  # noqa: C901
    r"""
    Produces an instrument definition file in XML format.

    Parameters
    ----------
    idf_xml_name: str
        Instrument geometry. Available options are:
        - 'rectangular detector' for a flat pixelated detector.
    request: ~pytest.request
        Parameters of the instrument.

    Returns
    -------
    str
    """
    #################
    # Below comes the functions in charge of creating specific instrument geometries
    #################
    def rectangular_detector(req_params):
        r"""
        Rectangular detector with Nx X Ny pixels

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
        """
        # use hidden attibutes to get data dimension, Nx and Ny can override this
        Nx, Ny, _ = _getDataDimensions(req_params)

        # get the parameters from the request object
        params = {'name': req_params.get('name', 'GenericSANS'),
                  'l1': -1. * abs(float(req_params.get('l1', -11.))),
                  'Nx': Nx,
                  'Ny': Ny,
                  'dx': float(req_params.get('dx', 1.)),
                  'dy': float(req_params.get('dy', 1.)),
                  'xcenter': float(req_params.get('xc', 0.)),
                  'ycenter': float(req_params.get('yc', 0.)),
                  'zcenter': float(req_params.get('zc', 5.))}
        params['dx_mm'] = params['dx'] * 1000.
        params['dy_mm'] = params['dy'] * 1000.

        # check that nothing is crazy
        assert (params['Nx'] > 1 and params['Nx'] < 300)
        assert (params['Ny'] >= 1 and params['Ny'] < 300)
        assert params['dx'] > 0.
        assert params['dy'] > 0.
        assert params['zcenter'] >= 0.

        # derived parameters
        params['half_dx'] = params['dx'] * .5
        params['half_dy'] = params['dy'] * .5
        params['xstart'] = -(params['Nx'] - 1) * params['half_dx']
        params['ystart'] = -(params['Ny'] - 1) * params['half_dy']

        template_xml = '''<?xml version="1.0" encoding="UTF-8"?>
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
            rot="180.0" axis-x="0" axis-y="1" axis-z="0">
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

    ###############
    # Below comes the calling to the specific instrument geometry constructor
    ###############
    # try to get the parent request, in case of sub-requests
    try:
        req_params = request.param
    except AttributeError:
        try:
            req_params = request._parent_request.param
        except AttributeError:
            req_params = dict()
    # Select the instrument geometry and apply the constructor
    constructor = {'rectangular detector': rectangular_detector}
    return constructor[idf_xml_name](req_params)


@pytest.fixture(scope='function')
def generic_IDF(request):
    return idf_xml_factory('rectangular detector', request)

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

    Nx, Ny, Naxis = _getDataDimensions(req_params)
    if y is not None:
        # force it to be a numpy array
        # this is a no-op if it is already the right type
        y = np.array(y)
        y = y.reshape((Nx, Ny, Naxis))
    else:
        y = np.zeros((Nx, Ny), dtype=float)
    y = y.ravel()

    if e is not None:
        e = np.array(e).ravel()
    else:
        e = np.sqrt(y)
        e[e == 0.] = 1.  # the default SANS likes
    if x is not None:
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


@pytest.fixture(scope='function')  # noqa: C901
def workspace_with_instrument(generic_IDF, request):
    r"""
    Workspace factory with a given instrument.

    The fixture is used by passing arguments to the `generic_IDF` fixture as a dictionary. For instance:

    .. :code:block:

       @pytest.mark.parametrize('workspace_with_instrument', [{'Nx': 3, 'Ny': 2}], indirect=True)

    Once inside the test function, the factory is called with the following optional arguments:
        output_workspace: Name of the workspace (default: random name prefixed with '__')
        axis_values : 1D array or list of the independent axis for the data. It will be copied across
            all spectra (default 0 for all spectra)
        intensities : ndarray or 2d/3d list of intensities for the instrument. Detector dimensions are inferred
            from the dimensionality. This will be linearized using `numpy.ravel`. (default: zeros of dimension Nx x Ny)
        uncertainties : ndarray or 2d/3d list of intensities for the instrument. This will be linearized using
            `numpy.ravel`. (default: sqrt(intensities), or one if intensity is zero)
        view: either 'array' or 'pixel'. In array-view the first index of the input data arrays travels each tube
            from top to bottom, and the second index travels across tubes. In pixel-view the first index travels
            across tubes and the second index travels each tube from bottom to top. (default 'pixel').
        axis_units : units for the independent axis (default 'wavelength')

    Example:
        ws2 = workspace_with_instrument(axis_values=[42.], intensities=[[1., 4.], [9., 16.], [25., 36.]])
    For more examples of use, look within testing class `test_fixtures.py::TestWorkspaceWithInstrument`

    Parameters
    ----------
    request: str
        A dictionary containing the following keys:
            name: Name of the instrument     (default: GenericSANS)
            Nx : number of columns (a.k.a number of tubes) (default 3)
            Ny : number of rows (a.k.a tube length)     (default 3)
            dx : width of a column in meters            (default 1)
            dy : height of a row in meters              (default 1)
            xc : distance of center along the x axis    (default 0)
            yc : distance of center along the y axis    (default 0)
            zc : distance of center along the z axis    (default 5)
            l1 : distance from source to sample       (default -11)

    Returns
    -------
    A function that can be used to generate multiple workspaces with the same instrument

    """
    try:
        instrument_params = request.param
    except AttributeError:
        try:
            instrument_params = request._parent_request.param
        except AttributeError:
            instrument_params = dict()

    # set up the default 'view' for the factory
    view = instrument_params.pop('view', 'pixel')

    workspace_inventory = list()  # holds created workspaces

    def factory(output_workspace=None, axis_units='wavelength',
                axis_values=None, intensities=None, uncertainties=None, view=view,
                number_x_pixels=instrument_params.get('Nx', 3), number_y_pixels=instrument_params.get('Ny', 3)):
        # Initialization of these options within the function signature results in the interpreter assigning a
        # function signature preserved through function call.
        if output_workspace is None:
            output_workspace = unique_workspace_dundername()

        if view not in ['array', 'pixel']:
            raise RuntimeError('Invalid value of view="{}". Must be "array" or "pixel"'.format(view))

        if intensities is not None:
            if isinstance(intensities, np.ndarray) is False:
                intensities = np.array(intensities)
            if view == 'array':
                if intensities.ndim == 2:
                    # reverse first index to increase tube ID along decreasing values on the X-axis
                    # reverse second index to increase pixel ID along each tube and along increasing values on Y-axis
                    intensities = np.transpose(intensities)[:, ::-1]
                elif intensities.ndim == 3:
                    intensities = np.transpose(intensities, axes=(1, 0, 2))[:, ::-1, :]
            number_x_pixels, number_y_pixels = intensities.shape[:2]
        else:
            intensities = np.zeros((number_x_pixels, number_y_pixels), dtype=float)
        intensities = intensities.ravel()

        if uncertainties is not None:
            uncertainties = np.array(uncertainties)
            if view == 'array':
                if uncertainties.ndim == 2:
                    uncertainties = uncertainties.transpose()[:, ::-1]
                elif uncertainties.ndim == 3:
                    uncertainties = np.transpose(uncertainties, axes=(1, 0, 2))[:, ::-1, :]
            uncertainties = uncertainties.ravel()
        else:
            uncertainties = np.sqrt(intensities)
            uncertainties[uncertainties == 0.] = 1.  # the default SANS likes

        if axis_values is not None:
            axis_values = np.array(axis_values)
        else:
            axis_values = np.zeros(1, dtype=float)

        n_pixels = number_x_pixels * number_y_pixels
        workspace = CreateWorkspace(DataX=axis_values, DataY=intensities, DataE=uncertainties, Nspec=n_pixels,
                                    UnitX=axis_units, OutputWorkspace=output_workspace)
        instrument_name = re.search(r'instrument name="([A-Za-z0-9_-]+)"', generic_IDF).groups()[0]
        LoadInstrument(Workspace=workspace, InstrumentXML=generic_IDF, RewriteSpectraMap=True,
                       InstrumentName=instrument_name)
        workspace_inventory.append(output_workspace)
        return workspace

    yield factory
    # Teardown
    for workspace_name in workspace_inventory:
        DeleteWorkspace(workspace_name)


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


def _assert_both_set_or_none(left, right, assert_func, err_msg):
    '''Either both argumentes are :py:obj:`None` or they are equal to each other'''
    if left is None and right is None:
        return
    if (left is not None) and (right is not None):
        assert_func(left, right, err_msg=err_msg)
    raise AssertionError('{}Either both or neither should be None (left={}, right={})'.format(err_msg, left, right))


def assert_wksp_equal(left, right, rtol=0, atol=0, err_msg=''):  # noqa: C901
    '''Generic method for checking equality of two data objects. This has some understanding of
    easily convertable types.'''
    id_left = getDataType(left)
    id_right = getDataType(right)

    # append colon to error message to make errors more readable
    if err_msg:
        err_msg += ': '

    # function pointer to make comparison code more flexible
    if rtol > 0 or atol > 0:
        assert_func = np.testing.assert_allclose
        kwargs = {'rtol': rtol, 'atol': atol}
    else:
        assert_func = np.testing.assert_equal
        kwargs = dict()

    # all of the comparison options - mixed modes first
    if id_left == DataType.WORKSPACE2D and id_right == DataType.IQ_MOD:
        units = left.getAxis(0).getUnit().caption()
        assert units == 'q', '{}: Found units="{}" rather than "q"'.format(err_msg, units)
        assert_func(left.extractX().ravel(), right.mod_q, err_msg=err_msg + 'mod_q', **kwargs)
        assert_func(left.extractY().ravel(), right.intensity, err_msg=err_msg + 'intensity', **kwargs)
        assert_func(left.extractE().ravel(), right.error, err_msg=err_msg + 'error', **kwargs)
    elif id_left == DataType.IQ_MOD and id_right == DataType.WORKSPACE2D:
        units = right.getAxis(0).getUnit().caption()
        assert units == 'q', '{}Found units="{}" rather than "q"'.format(err_msg, units)
        assert_func(left.mod_q, right.extractX().ravel(), err_msg=err_msg + 'mod_q', **kwargs)
        assert_func(left.intensity, right.extractY().ravel(), err_msg=err_msg + 'intensity', **kwargs)
        assert_func(left.error, right.extractE().ravel(), err_msg=err_msg + 'error', **kwargs)
    elif id_left == id_right:  # compare things that are the same type
        if id_left == DataType.WORKSPACE2D:
            # let mantid do all the work
            if atol > 0:
                cmp, messages = CompareWorkspaces(Workspace1=str(left), Workspace2=str(right), Tolerance=atol)
            else:
                cmp, messages = CompareWorkspaces(Workspace1=str(left), Workspace2=str(right),
                                                  Tolerance=rtol, ToleranceRelErr=True)
            messages = [row['Message'] for row in messages]
            assert cmp, err_msg + '; '.join(messages)
        else:
            # all the other data objects share some attributes
            assert_func(left.intensity, right.intensity, err_msg=err_msg + 'intensity', **kwargs)
            assert_func(left.error, right.error, err_msg=err_msg + 'error', **kwargs)
            _assert_both_set_or_none(left.wavelength, right.wavelength, assert_func, err_msg + 'wavelength')
            if id_left == DataType.IQ_MOD:
                assert_func(left.mod_q, right.mod_q, err_msg=err_msg + 'mod_q', **kwargs)
                _assert_both_set_or_none(left.delta_mod_q, right.delta_mod_q, assert_func, err_msg + 'delta_mod_q')
            elif id_left == DataType.IQ_AZIMUTHAL:
                assert_func(left.qx, right.qx, err_msg=err_msg + 'qx', **kwargs)
                assert_func(left.qy, right.qy, err_msg=err_msg + 'qy', **kwargs)
                _assert_both_set_or_none(left.delta_qx, right.delta_qx, assert_func, err_msg + 'delta_qx')
                _assert_both_set_or_none(left.delta_qy, right.delta_qy, assert_func, err_msg + 'delta_qy')
            elif id_left == DataType.IQ_CRYSTAL:
                assert_func(left.qx, right.qx, err_msg=err_msg + 'qx', **kwargs)
                assert_func(left.qy, right.qy, err_msg=err_msg + 'qy', **kwargs)
                assert_func(left.qz, right.qz, err_msg=err_msg + 'qz', **kwargs)
                _assert_both_set_or_none(left.delta_qx, right.delta_qx, assert_func, err_msg + 'delta_qx')
                _assert_both_set_or_none(left.delta_qy, right.delta_qy, assert_func, err_msg + 'delta_qy')
                _assert_both_set_or_none(left.delta_qz, right.delta_qz, assert_func, err_msg + 'delta_qz')
            else:
                raise NotImplementedError('Do not know how to compare {} objects'.format(id_left))
    else:
        raise NotImplementedError('Do not know how to compare {} and {}'.format(id_left, id_right))


def depth(L):
    '''
    calculating the depth of the object
    :param L: the given object (float, integer, list, list of list)
    :return: integer value of depth
    '''
    return (isinstance(L, list) and max(map(depth, L)) + 1)


def pixel_block(pixel_number, radius, height):
    '''
    generating the cylindrical shape of the detector
    :param pixel_number:integer of the number of pixels
    :param radius: float of the radius of cylinder detector in meter
    :param height: float of the height of the cylinder detector in meter
    :return: the .xml format generating a cylindrical shape of the detector from height and radius
    '''
    return """
    <type is="detector" name="pixel_{pixel_number}">
        <cylinder id="cyl-approx">
            <centre-of-bottom-base x="0" y="-{half_y}" z="0"/>
            <axis x="0.0" y="1.0" z="0.0"/>
            <radius val="{radius}"/>
            <height val="{height}"/>
        </cylinder>
        <algebra val="cyl-approx"/>
    </type>
""".format(pixel_number=pixel_number, half_y=height*0.5, radius=radius, height=height)


def pixel_location(pixel_number, xcenter, ycenter, zcenter):
    '''
    specifying the pixel location in the space
    :param pixel_number: integer of the number of pixels
    :param xcenter:distance of center along the x axis
    :param ycenter:distance of center along the y axis
    :param zcenter:distance of center along the z axis
    :return: the .xml format of specifying the pixel location in the space
    '''
    return """
        <component type="pixel_{pixel_number}" name="pixel_{pixel_number}">
            <location y="{ycenter}" x="{xcenter}" z="{zcenter}"/>
        </component>
""".format(pixel_number=pixel_number, xcenter=xcenter, ycenter=ycenter, zcenter=zcenter)


@pytest.fixture(scope='function')
def arbitrary_assembly_IDF(request):
    '''
    generate a test IDF with a cylindrical detector pixels

    Parameters
    ----------

    request is a dictionary containing the following keys:

        name: Name of the instrument     (default: GenericSANS)
        radius : list of radius  in meters            (default 1)
        height : list of height in meters              (default 1)
        pixel_center : list of list of distance of center along the x axis,
                      along the y axis and along the z axis
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

    # get the parameters from the request object
    params = {'name': req_params.get('name', 'GenericSANS'),
              'l1': -1. * abs(float(req_params.get('l1', -11.))),
              'radius': req_params.get('radius', 0.0025),
              'height': req_params.get('height', 0.005),
              'pixel_centers': req_params.get('pixel_centers'),
              }

    number_pixels = len(params['pixel_centers'])

    # check that nothing is crazy
    if (depth(params['radius'])) == 0:
        params['radius'] = [params['radius']]*number_pixels
    if (depth(params['height'])) == 0:
        params['height'] = [params['height']] * number_pixels
    assert (len(params['radius']) > 1 and len(params['radius']) < 300)
    assert (len({len(params['radius']), len(params['height']),
                 len(params['pixel_centers'])}) == 1)
    assert params['radius'] > [0]*number_pixels
    assert params['height'] > [0]*number_pixels
    assert list(list(zip(*params['pixel_centers']))[-1]) >= [0]*number_pixels  # zcenter has to be positive

    pixel_blocks = [pixel_block(i, params['radius'][i],
                                params['height'][i]) for i in range(number_pixels)]
    pixel_blocks = '\n'.join(pixel_blocks)

    pixel_locations = [pixel_location(i, *params['pixel_centers'][i])
                       for i in range(number_pixels)]

    pixel_locations = '\n'.join(pixel_locations)

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


    <!-- Pixel for Detectors-->
{pixel_blocks}

    <type name="arbitrary_assembly">
{pixel_locations}
    </type>

     <!-- Arbitrary Assembly of Cylindrical Detector-->
    <component type="arbitrary_assembly" name="detector1" idlist="pixel_ids">
        <location/>
    </component>

  <!---->
  <!--LIST OF PIXEL IDs in DETECTOR-->
  <!---->
  <idlist idname="pixel_ids">
    <id end="{max_pixels_id}" start="0"/>
    </idlist>
</instrument>'''

    # return the completed template
    return template_xml.format(name=params['name'], l1=params['l1'],
                               pixel_blocks=pixel_blocks, pixel_locations=pixel_locations,
                               max_pixels_id=number_pixels-1)
