from mantid.api import AnalysisDataService, FileFinder
from mantid.kernel import ConfigService
from drtsans.instruments import instrument_label, extract_run_number
import os
import stat
import pathlib

__all__ = ['abspath', 'exists', 'registered_workspace', 'allow_overwrite']


def allow_overwrite(folder):
    r"""
    Changes permissions for all the files and folders in the path
    to allow write by anyone. It is not recursive. It will do that
    only for the files/folders the user has permissions to do that

    Parameters
    ----------
    path: str
        string to the folder in which the file permissions are to be changed

    Returns
    -------
    None
    """
    for path in pathlib.Path(folder).glob('*'):
        permissions = path.stat().st_mode
        try:
            path.chmod(permissions | stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
        except PermissionError:
            pass


def abspath(path, instrument='', ipts=''):
    r"""
    Returns an absolute path

    In addition to fully supporting what os.path.abspath does,
    this also supports path strings in such as ``EQSANS_106026`` and
    ``EQSANS106026``. It will search your data search path and the
    data archive using ONCat.

    This uses mantid.api.FileFinder.
    """
    # don't use network for first check
    if os.path.exists(path):
        return os.path.abspath(path)

    # guess the path from existing information
    if ipts:
        runnumber = path.lower()
        if not instrument:
            try:
                instrument = instrument_label(path)
                runnumber = extract_run_number(path)
            except RuntimeError:
                pass  # failed to extract instrument/runnumber

        if instrument:  # only try if instrument is known
            try:
                runnumber = int(runnumber)  # make sure it doesn't have extra characters
                instrument = ConfigService.getInstrument(instrument)  # to object
                facility = str(instrument.facility())
                instrument = str(instrument)  # convert back to short name
                option = f'/{facility}/{instrument}/IPTS-{ipts}/nexus/{instrument}_{runnumber}.nxs.h5'
                print('Seeing if "{}" exists'.format(option))
                if os.path.exists(option):
                    return option
            except RuntimeError:
                pass  # facility not found
            except ValueError:
                pass  # could not convert run number to an integer

            # prepend the instrument name if it isn't already there
            if not path.upper().startswith(instrument.upper()):
                path = '{}{}'.format(instrument.upper(), path)

    # get a full path from `datasearch.directories`
    option = FileFinder.getFullPath(path)
    if option and os.path.exists(option):
        return option

    # get all of the options from FileFinder and convert them to an absolute
    # path in case any weren't already
    try:
        options = [os.path.abspath(item) for item in FileFinder.findRuns(path)]
    except RuntimeError:
        options = []  # marks things as broken
    if not options:  # empty result
        raise RuntimeError('Failed to find location of file from hint '
                           '"{}"'.format(path))
    for option in options:
        if os.path.exists(option):
            return option

    raise RuntimeError('None of the locations suggested by ONCat contain '
                       'existing files for "{}"'.format(path))


def exists(path):
    r"""
    Test whether a path exists.  Returns False for broken symbolic links

    In addition to fully supporting what os.path.exists does,
    this also supports path strings in such as ``EQSANS_106026`` and
    ``EQSANS106026``. It will search your data search path and the
    data archive using ONCat.

    This uses mantid.api.FileFinder.
    """
    # quickest way is to assume it is a regular file
    if os.path.exists(path):
        return True
    # check in `datasearch.directories`
    if os.path.exists(FileFinder.getFullPath(path)):
        return True
    # check via locations provided by ONCat
    try:
        for option in FileFinder.findRuns(path):
            if os.path.exists(option):
                return True
    except RuntimeError:
        return False  # no suggestions found

    return False


def registered_workspace(source):
    r"""
    Find out if the source is a workspace registered in the Analysis Data
    Service.

    Parameters
    ----------
    source: str, Workspace

    Returns
    -------
    bool
    """
    return AnalysisDataService.doesExist(str(source))
