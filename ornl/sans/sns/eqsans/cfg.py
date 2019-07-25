"""
Reader for EQSANS configuration files in the old format
"""

from __future__ import (absolute_import, division, print_function)

import os
import re
import numpy as np
from contextlib import contextmanager

cfg_dir = '/SNS/EQSANS/shared/instrument_configuration'


def closest_config(run, config_dir=cfg_dir):
    """
    Configuration file for a given run number.
    The appropriate configuration file is the one marked with a run number
    that is closest and smaller (or equal) than the input run number

    Parameters
    ----------
    config_dir: str
        Directory containing configuration files

    Returns
    -------
    str
        Absolute path to configuration file
    """
    pattern = re.compile(r'eqsans_configuration\.(\d+)')
    reference_runs = list()
    for root, dirs, files in os.walk(config_dir):
        for file_name in files:
            match = pattern.search(file_name)
            if match is not None:
                reference_runs.append(int(match.groups()[0]))
    reference_runs.sort()
    reference_runs = np.asarray(reference_runs)
    maximum_index_below_run = np.where(run >= reference_runs)[0][-1]
    reference_run = reference_runs[maximum_index_below_run]
    return os.path.join(config_dir,
                        'eqsans_configuration.{}'.format(reference_run))


@contextmanager
def open_source(source, config_dir=cfg_dir):
    """
    Find the configuration file appropriate to the input source info

    Parameters
    ----------
    source: str
        One of the following: (1) absolute path or just filename to a
        configuration file; (2) run-number
    config_dir: str
        Directory containing configuration files

    Yields
    ------
    file handle:
        Handle to the configuration file
    """
    src = str(source)
    if os.path.isfile(src):
        file_handle = open(src)
    else:
        file = os.path.join(config_dir, src)
        if os.path.isfile(file):
            file_handle = open(file)
        else:
            run = int(source)
            file_handle = open(closest_config(run, config_dir=config_dir))
    try:
        yield file_handle
    finally:
        file_handle.close()


class CfgItemValue(object):
    """
    Entry value in EQSANS configuration file

    Parameters
    ----------
    data: string, or list of strings
        value of the entry
    off: bool
        True if the entry was commented in the file
    note: str
        Description of the entry
    """

    def __init__(self, name='', data='', off=True, note=''):
        self.data = data
        self.off = off
        self.note = note

    def __repr__(self):
        return 'CfgItemValue(data="{data}", off={off}, note="{note}")'.\
            format(**self.__dict__)

    def __eq__(self, other):
        """Discard note explanatory when comparing two value items"""
        return self.data == other.data and self.off == other.off


class Cfg(object):
    """
    Read EQSANS configuration files
    """

    @staticmethod
    def load(source, config_dir=cfg_dir):
        """
        Load the configuration file appropriate to the input source info

        Parameters
        ----------
        source: str
            One of the following: (1) absolute path or just filename to a
            configuration file; (2) run-number
        config_dir: str
            Directory containing configuration files

        Returns
        -------
        dict
            A dictionary with CfgItemValue objects as values
        """
        cfg = dict()
        with open_source(source, config_dir=config_dir) as f:
            for line in f.readlines():
                if '=' not in line:
                    continue  # this line contains no valid entries
                key, val = [x.strip() for x in line.split('=')]
                commented = False
                description = ''
                if '#' in val:
                    val, description = [x.strip() for x in val.split('#')]
                if '#' in key:
                    commented = True
                    key = key.split('#')[-1].strip()
                if key in cfg:
                    if commented is True:
                        continue
                    old_val = cfg[key].data
                    if isinstance(old_val, list):
                        cfg[key].data.append(val)
                    else:
                        cfg[key].data = [old_val, val]
                    if description != '':
                        cfg[key].help = description
                else:
                    item = CfgItemValue(data=val, off=commented,
                                        note=description)
                    cfg[key] = item
        return cfg

    def __init__(self, source=None, config_dir=cfg_dir):
        self._cfg = dict() if source is None \
            else Cfg.load(source, config_dir=config_dir)

    def __getitem__(self, item):
        return self._cfg[item]

    def __setitem__(self, key, value):
        self._cfg[key] = value

    def as_dict(self):
        """Return all data as a dictionary of key data values"""
        return {k: v.data for (k, v) in self._cfg.items()}

    def __repr__(self):
        fmt = '"{}" : {}'
        return '\n'.join(fmt.format(k, v) for (k, v) in self._cfg.items())

    def logit(self, key, workspace, name=None, replace=False):
        """

        Parameters
        ----------
        key: str
            Key associated to a specific configuration entry.
        workspace: mantid.MatrixWorkspace
            Save the property in the logs of this workspace
        name: str
            Alternative log name to key
        replace: bool
            Overwrite existing log entry

        Raises
        ------
        ValueError
            If the log entry exists and replace is set to False
        """
        log_name = key if name is None else name
        run = workspace.getRun()
        if replace is False and run.hasProperty(log_name):
            msg = 'Property {} already exists. Set keyword "replace"' \
                  ' to True if you wish to replace the existing property' \
                  ' with the new value.'
            raise ValueError(msg.format(log_name))
        run.addProperty(log_name, self[key].data, replace=replace)
