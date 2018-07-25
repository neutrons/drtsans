#!/usr/bin/env python

from __future__ import print_function

'''
Run as simple test

PYTHONPATH=. pytest -v -s tests/test_eqsansload.py

-s : Shows the std out

'''
import sys
import os

from dotenv import load_dotenv
load_dotenv()

def test_get_config_file():
    from ornl.sans.sns.eqsans.parameters import _get_config_file
    assert _get_config_file(71820) == '/SNS/EQSANS/shared/instrument_configuration/eqsans_configuration.71820'
    assert _get_config_file(71821) == '/SNS/EQSANS/shared/instrument_configuration/eqsans_configuration.71820'
    assert _get_config_file(72001) == '/SNS/EQSANS/shared/instrument_configuration/eqsans_configuration.71820'


def test_get_parameters():
    from ornl.sans.sns.eqsans.parameters import get_parameters
    params = get_parameters(68200)
    assert params['detector pixel sizes'] == '5.5, 5.5'
    assert params['rectangular mask'].split('\n')[0] == '0 0;7 255'


def test_EQSANSLoad():
    '''
    EQSANSLoad workflow algorithm as called by Mantid
    '''
    
    sys.path.append(os.getenv("MANTID_PATH"))
    from mantid.simpleapi import EQSANSLoad
    from mantid.kernel import PropertyManagerDataService, PropertyManager

    pm = PropertyManager()
    PropertyManagerDataService.addOrReplace("test_pm", pm)
    out = EQSANSLoad(
        Filename=os.path.join(os.getenv('DATA_DIRECTORY'), 'eqsans', 'EQSANS_68200_event.nxs'),
        # UseDirectBeamMethod=True,
        # BeamRadius=3,
        ReductionProperties='test_pm',
    )

    x = float(pm.getPropertyValue('LatestBeamCenterX'))
    y = float(pm.getPropertyValue('LatestBeamCenterY'))

    assert x == 95.5
    assert y == 127.5

    print(out.OutputMessage)
