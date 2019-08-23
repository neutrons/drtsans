from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.mark.parametrize('generic_instrument', [{'Nx': 2, 'Ny': 2}], indirect=True)
def test_dark_current_subtract_dark_current(generic_instrument):
    ws = generic_instrument
    print(ws.extractX())
    assert False


if __name__ == '__main__':
    pytest.main()
