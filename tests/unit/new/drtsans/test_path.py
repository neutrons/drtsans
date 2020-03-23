import pytest
import pathlib
import stat
from mantid.simpleapi import CreateWorkspace
from drtsans.path import abspath, exists, registered_workspace, allow_overwrite
from drtsans.settings import amend_config, unique_workspace_dundername as uwd
from os.path import exists as os_exists
from tempfile import gettempdir, NamedTemporaryFile

# arbitrarily selected IPTS to see if the mount is in place
HAVE_EQSANS_MOUNT = os_exists('/SNS/EQSANS/IPTS-23732/')

SEARCH_ON = {}
SEARCH_OFF = {'datasearch.searcharchive': 'off'}

IPTS_23732 = '/SNS/EQSANS/IPTS-23732/nexus/'


@pytest.mark.skipif(not HAVE_EQSANS_MOUNT,
                    reason='Do not have /SNS/EQSANS properly '
                    'mounted on this system')
@pytest.mark.parametrize('hint, fullpath',
                         [('EQSANS_106026',
                           IPTS_23732 + 'EQSANS_106026.nxs.h5'),
                          ('EQSANS106027',
                           IPTS_23732 + 'EQSANS_106027.nxs.h5'),
                          ('EQSANS_88974.nxs.h5',
                           'DATADIR/EQSANS_88974.nxs.h5')],
                         ids=('EQSANS_106026', 'EQSANS_106026',
                              'EQSANS_88974'))
def test_abspath_with_archivesearch(hint, fullpath, reference_dir):
    # set the data directory in the result using the test fixture
    fullpath = fullpath.replace('DATADIR', reference_dir.new.eqsans)

    with amend_config(SEARCH_ON, data_dir=reference_dir.new.eqsans):
        assert abspath(hint) == fullpath


@pytest.mark.parametrize('hint',
                         ['randomname', 'EQSANS_106026', 'EQSANS_106026'],
                         ids=('randomname', 'EQSANS_106026', 'EQSANS_106026'))
def test_abspath_without_archivesearch(hint):
    with amend_config(SEARCH_OFF):
        with pytest.raises(RuntimeError):
            found = abspath(hint)
            assert False, 'found "{}" at "{}"'.format(hint, found)


@pytest.mark.skipif(not HAVE_EQSANS_MOUNT,
                    reason='Do not have /SNS/EQSANS properly '
                    'mounted on this system')
@pytest.mark.parametrize('hint, found',
                         [('EQSANS_106026', True),
                          ('EQSANS106027', True),
                          ('EQSANS_88974.nxs.h5', True)])
def test_exists_with_archivesearch(hint, found, reference_dir):
    with amend_config(SEARCH_ON, data_dir=reference_dir.new.eqsans):
        assert exists(hint) == found  # allows verifying against True and False


@pytest.mark.parametrize('hint, found',
                         [('EQSANS_106026', False),
                          ('EQSANS106027', False),
                          ('EQSANS_88974.nxs.h5', True)])
def test_exists_without_archivesearch(hint, found, reference_dir):
    with amend_config(SEARCH_OFF, data_dir=reference_dir.new.eqsans):
        assert exists(hint) == found  # allows verifying against True and False


def test_registered_workspace():
    w_name = uwd()
    assert registered_workspace(w_name) is False
    w = CreateWorkspace(DataX=[1], Datay=[1], OutputWorkspace=w_name)
    assert registered_workspace(w_name) is True
    assert registered_workspace(w) is True


def test_allow_overwrite():
    tmpdir = gettempdir()
    # create an empty file
    tmpfile = NamedTemporaryFile(dir=tmpdir, delete=False)
    tmpfile.close()
    # check if others write permission is false
    path = pathlib.Path(tmpfile.name)
    assert not bool(path.stat().st_mode & stat.S_IWOTH)
    allow_overwrite(tmpdir)
    # check permissions
    assert bool(path.stat().st_mode & stat.S_IWUSR)
    assert bool(path.stat().st_mode & stat.S_IWGRP)
    assert bool(path.stat().st_mode & stat.S_IWOTH)
    # delete file
    path.unlink()


if __name__ == '__main__':
    pytest.main([__file__])
