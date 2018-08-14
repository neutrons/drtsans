from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.mark.skip(reason="only for debugging")
def test_eqsans_w(eqsans_w):
    pass


#@pytest.mark.skip(reason="only for debugging")
def test_porasil_slice1m(porasil_slice1m):
    for k in porasil_slice1m.w.keys():
        assert w._w[k].name() == '_'+k
        assert w[k].name() == k


if __name__ == '__main__':
    pytest.main()
