from __future__ import (absolute_import, division, print_function)

import pytest
import os
from ornl.sans.sns.eqsans import cfg


def test_setitem():
    c = cfg.Cfg()
    value = cfg.CfgItemValue(data=42, off=False, note='meaning of universe')
    c['k'] = value  # set value
    assert c['k'] == value  # get value


def test_closest_config():
    name = 'eqsans_configuration.92474'
    assert cfg.closest_config(97711) == os.path.join(cfg.cfg_dir, name)


def test_open_source():
    name = 'eqsans_configuration.92474'
    full_name = os.path.join(cfg.cfg_dir, name)
    with cfg.open_source(full_name) as f:
        assert f.name == full_name
    with cfg.open_source(name) as f:
        assert f.name == full_name
    with cfg.open_source(97711) as f:
        assert f.name == full_name


def test_load():
    c = cfg.Cfg(source=97711)
    value = cfg.CfgItemValue(data='500 2000', off=False)
    assert c['TOF edge discard'] == value


if __name__ == '__main__':
    pytest.main()
