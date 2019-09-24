from __future__ import (absolute_import, division, print_function)

import pytest


@pytest.mark.parametrize('generic_workspace', [{'Nx': 2, 'Ny': 2}], indirect=True)
def test_dark_current_subtract_dark_current(generic_workspace):
    ws = generic_workspace
    print(ws.extractX())
    assert True


if __name__ == '__main__':
    pytest.main()
