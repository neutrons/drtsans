from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.mark.skip(reason="only for debugging")
def test_eqsans_w(eqsans_w):
    pass


@pytest.mark.skip(reason="only for debugging")
def test_porasil_slice1m(porasil_slice1m):
    w = porasil_slice1m.w
    for k in w.keys():
        assert w._w[k].name() == '_porasil_slice1m_' + k
        assert w[k].name() == 'porasil_slice1m_' + k


@pytest.mark.skip(reason="only for debugging")
def test_frame_skipper(frame_skipper):
    w = frame_skipper.w
    for k in w.keys():
        assert w._w[k].name() == '_frame_skipper_' + k
        assert w[k].name() == 'frame_skipper_' + k


if __name__ == '__main__':
    pytest.main()
