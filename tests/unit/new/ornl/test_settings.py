from __future__ import (absolute_import, division, print_function)

import pytest
from ornl.settings import (namedtuplefy, amend_config,
                           unique_workspace_dundername,
                           optional_output_workspace)
from mantid.simpleapi import CreateWorkspace
from mantid.kernel import ConfigService


@pytest.mark.offline
def test_namedtuplefy():

    @namedtuplefy
    def foo(x):
        return dict(foo=x)

    @namedtuplefy
    def goo(x):
        return dict(goo=x)

    y1 = foo(42)
    z1 = goo(24)
    y2 = foo(41)
    z2 = goo(21)

    assert type(y1) == type(y2)
    assert type(z1) == type(z2)
    assert type(y1) != type(z1)
    assert type(y2) != type(z2)


@pytest.mark.offline
def test_amend_config():
    config = ConfigService.Instance()
    old_instrument = config['instrumentName']
    with amend_config({'instrumentName': '42'}):
        assert config['instrumentName'] == '42'
    assert config['instrumentName'] == old_instrument


@pytest.mark.offline
def test_offline():
    print('this tests runs when offline')


def test_optional_output_workspace():

    @optional_output_workspace
    def foo():
        r"""Creates a silly workspace with some random name"""
        return CreateWorkspace(DataX=[42, ], DataY=[42, ],
                               OutputWorkspace=unique_workspace_dundername())
    ws = foo()
    assert ws.name() == 'ws'
    ws = foo(output_workspace='sw')
    assert ws.name() == 'sw'

    @optional_output_workspace
    def foo(output_workspace=None):
        r"""Creates a silly workspace with some random name"""
        return CreateWorkspace(DataX=[42, ], DataY=[42, ],
                               OutputWorkspace=output_workspace)
    ws = foo()
    assert ws.name() == 'ws'
    ws = foo(output_workspace='sw')
    assert ws.name() == 'sw'

    @optional_output_workspace
    def foo(output_workspace='foo_name'):
        r"""Creates a silly workspace with some random name"""
        return CreateWorkspace(DataX=[42, ], DataY=[42, ],
                               OutputWorkspace=output_workspace)
    ws = foo()
    assert ws.name() == 'foo_name'
    ws = foo(output_workspace='sw')
    assert ws.name() == 'sw'

    @optional_output_workspace
    def foo(iws):
        return 2.0 * iws
    ws = CreateWorkspace(DataX=[42, ], DataY=[42, ],
                         OutputWorkspace='meaning_of_the_universe')
    other_ws = foo(ws, output_workspace='skeletor')
    assert other_ws.name() == 'skeletor'
    other_ws = foo(ws)
    assert other_ws.name() == 'meaning_of_the_universe'

    # Corner case: returned workspace is the input workspace
    @optional_output_workspace
    def foo(iws):
        return iws
    ws = CreateWorkspace(DataX=[42, ], DataY=[42, ],
                         OutputWorkspace='meaning_of_the_universe')
    other_ws = foo(ws, output_workspace='skeletor')
    assert other_ws.name() == 'skeletor'  # a simple Rename
    other_ws = foo(ws)
    assert other_ws.name() == 'skeletor'  # nothing is done


if __name__ == '__main__':
    pytest.main()
