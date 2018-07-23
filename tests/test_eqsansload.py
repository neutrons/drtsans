#!/usr/bin/env python

from __future__ import print_function
'''
Run as simple test

PYTHONPATH=. pytest tests/test_eqsansload.py

'''

from configparser import RawConfigParser
from itertools import chain

from ornl.settings import MultiOrderedDict


def test_read_configuration():
    conf_file = "/SNS/EQSANS/shared/instrument_configuration/eqsans_configuration.92474"

    parser = RawConfigParser(
        dict_type=MultiOrderedDict,
        strict=False,
        inline_comment_prefixes="#")
    with open(conf_file, 'r') as f:
        f = chain(("[DEFAULT]",), f)  # This line does the trick.
        parser.read_file(f)
    
    assert parser['DEFAULT']['rectangular mask'].split("\n")[0] == "0 0;0 255"
