from mantid.api import FileFinder
from os import path as os_path

__all__ = ['abspath', 'exists']


def abspath(path):
    '''
    Returns an absolute path

    This uses :ref:`mantid.api.FileFinder`.
    '''
    # don't use network for first check
    if os_path.exists(path):
        return os_path.abspath(path)

    # get a full path from `datasearch.directories`
    option = FileFinder.getFullPath(path)
    if option and os_path.exists(option):
        return option

    # get all of the options from FileFinder and convert them to an absolute
    # path in case any weren't already
    try:
        options = [os_path.abspath(item) for item in FileFinder.findRuns(path)]
    except RuntimeError:
        options = []  # marks things as broken
    if not options:  # empty result
        raise RuntimeError('Failed to find location of file from hint '
                           '"{}"'.format(path))
    for option in options:
        if os_path.exists(option):
            return option

    raise RuntimeError('None of the locations suggested by ONCat contain '
                       'existing files for "{}"'.format(path))


def exists(path):
    '''
    Test whether a path exists.  Returns False for broken symbolic links

    This uses :ref:`mantid.api.FileFinder`.
    '''
    # quickest way is to assume it is a regular file
    if os_path.exists(path):
        return True
    # check in `datasearch.directories`
    if os_path.exists(FileFinder.getFullPath(path)):
        return True
    # check via locations provided by ONCat
    try:
        for option in FileFinder.findRuns(path):
            if os_path.exists(option):
                return True
    except RuntimeError:
        return False  # no suggestions found

    return False
