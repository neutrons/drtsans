#!/usr/bin/env python

from __future__ import print_function

'''
Run as simple test

PYTHONPATH=. pytest -v -s tests/test_eqsansload.py

-s : Shows the std out

'''


def test_get_config_file():
    from ornl.sans.sns.eqsans.parameters import _get_config_file
    fn = '/SNS/EQSANS/shared/instrument_configuration/eqsans_configuration.{}'
    assert _get_config_file(71820) == fn.format(71820)
    assert _get_config_file(71821) == fn.format(71820)
    assert _get_config_file(72001) == fn.format(71820)


def test_get_parameters():
    from ornl.sans.sns.eqsans.parameters import get_parameters
    params = get_parameters(68200)
    assert params['detector pixel sizes'] == '5.5, 5.5'
    assert params['rectangular mask'].split('\n')[0] == '0 0;7 255'


def test_EQSANSLoad(eqsans_f):
    """
    EQSANSLoad workflow algorithm as called by Mantid
    """

    from mantid.simpleapi import EQSANSLoad
    from mantid.kernel import PropertyManagerDataService, PropertyManager

    pm = PropertyManager()
    PropertyManagerDataService.addOrReplace("test_pm", pm)
    out = EQSANSLoad(
        Filename=eqsans_f["darkcurrent"],
        # UseDirectBeamMethod=True,
        # BeamRadius=3,
        ReductionProperties='test_pm',
    )

    x = float(pm.getPropertyValue('LatestBeamCenterX'))
    y = float(pm.getPropertyValue('LatestBeamCenterY'))

    assert x == 95.5
    assert y == 127.5

    print(out.OutputMessage)
