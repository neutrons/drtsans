from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.mark.skip(reason="only for debugging")
def test_eqsans_w(eqsans_w):
    pass


@pytest.mark.skip(reason="only for debugging")
def test_porasil_slice1m(porasil_slice1m):
    assert porasil_slice1m()['w'] == dict()
    assert set(porasil_slice1m('dc')['w'].keys()) == set(('dc',))
    assert set(porasil_slice1m(('dc', 's'))['w'].keys()) == set(('dc', 's'))
    info_dict = porasil_slice1m('all')
    assert set(info_dict['w'].keys()) == set(info_dict['f'].keys())


if __name__ == '__main__':
    pytest.main()
