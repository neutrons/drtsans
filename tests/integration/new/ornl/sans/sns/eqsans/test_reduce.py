from __future__ import (absolute_import, division, print_function)

import pytest
from os.path import join as pj
from mantid.simpleapi import LoadNexus, SumSpectra, CompareWorkspaces

from ornl.sans.sns.eqsans import reduce
from ornl.settings import uwn


def test_load_w(refd):
    _w0 = reduce.load_w('EQSANS_92353', uwn(), ltc=500, htc=2000, dw=0.1)
    _w1 = SumSpectra(_w0)
    fn = pj(refd.new.eqsans, 'test_reduce', 'compare', 'ref_load_w.nxs')
    _w2 = LoadNexus(fn, OutputWorkspace=uwn())
    assert CompareWorkspaces(_w1, _w2)


if __name__ == '__main__':
    pytest.main()
