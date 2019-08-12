from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pj
from ornl.sans.sns.eqsans import cfg


@pytest.mark.offline
def test_setitem():
    c = cfg.Cfg()
    value = cfg.CfgItemValue(data=42, off=False, note='meaning of universe')
    c['k'] = value  # set value
    assert c['k'] == value  # get value


def test_closest_config(refd):
    config_dir = pj(refd.new.eqsans, 'instrument_configuration')
    name = pj(config_dir, 'eqsans_configuration.92474')
    assert cfg.closest_config(97711, config_dir=config_dir) == name


def test_open_source(refd):
    config_dir = pj(refd.new.eqsans, 'instrument_configuration')
    name = 'eqsans_configuration.92474'
    full_name = pj(config_dir, name)
    with cfg.open_source(full_name) as f:
        assert f.name == full_name
    with cfg.open_source(name, config_dir=config_dir) as f:
        assert f.name == full_name
    with cfg.open_source(97711, config_dir=config_dir) as f:
        assert f.name == full_name


def test_load(refd):
    config_dir = pj(refd.new.eqsans, 'instrument_configuration')
    c = cfg.Cfg(source=97711, config_dir=config_dir)
    value = cfg.CfgItemValue(data='500 2000', off=False)
    assert c['TOF edge discard'] == value
    assert isinstance(c['Rectangular Mask'], cfg.CfgItemRectangularMask)
    d = c.as_dict()
    assert d['Rectangular Mask'] == c['Rectangular Mask'].detectors


def test_rec_mask():
    # second rectangle is included in the first
    c = cfg.CfgItemRectangularMask(data=['0, 0; 1, 255', '1 ,250; 1, 255'])
    px = c.pixels
    assert (1, 255) in px


def test_mask_mixin():
    c = cfg.CfgItemRectangularMask(data=['1, 0; 2, 255', '1 250 1 255'])
    dets = c.detectors
    assert len(dets) == 2 * 256
    assert dets[-1] == 3 * 256 - 1
    assert c.value == dets


if __name__ == '__main__':
    pytest.main()
