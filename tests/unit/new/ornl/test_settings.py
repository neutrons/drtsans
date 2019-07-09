from __future__ import (absolute_import, division, print_function)

import pytest
from ornl.settings import (namedtuplefy, amend_config,
                           optional_output_workspace)
from mantid.simpleapi import mtd, CreateWorkspace, RenameWorkspace
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
    def foo(**kwargs):
        r"""Creates a silly workspace with some random name"""
        return CreateWorkspace(DataX=[42, ], DataY=[42, ],
                               OutputWorkspace=kwargs['output_workspace'])
    ws = foo()
    assert ws.name().startswith('__')
    ws = foo(output_workspace='sw')
    assert ws.name().startswith('sw')

    @optional_output_workspace
    def foo(output_workspace=None):
        r"""Creates a silly workspace with some random name"""
        return CreateWorkspace(DataX=[42, ], DataY=[42, ],
                               OutputWorkspace=output_workspace)
    ws = foo()
    assert ws.name().startswith('__')
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
    def foo(iws, output_workspace='meaning_of_the_universe', **kwargs):
        tmp = 2.0 * iws
        RenameWorkspace(InputWorkspace=tmp, OutputWorkspace=output_workspace)
        return mtd[output_workspace]

    ws = CreateWorkspace(DataX=[42, ], DataY=[42, ],
                         OutputWorkspace='meaning_of_the_universe')
    other_ws = foo(ws, output_workspace='skeletor')
    assert other_ws.name() == 'skeletor'
    other_ws = foo(ws)
    assert other_ws.name() == 'meaning_of_the_universe'

    @optional_output_workspace
    def foo(input_workspace, **kwargs):
        input_workspace = str(input_workspace)
        output_workspace = str(kwargs['output_workspace'])
        if input_workspace == output_workspace:
            return mtd[input_workspace]
        else:
            RenameWorkspace(InputWorkspace=input_workspace,
                            OutputWorkspace=output_workspace)
            return mtd[output_workspace]

    # test in-place as the default
    ws = CreateWorkspace(DataX=[42, ], DataY=[42, ],
                         OutputWorkspace='meaning_of_the_universe')
    other = foo(ws)
    assert other.name() == ws.name()
    other = foo(ws, output_workspace='skeletor')
    assert other.name() == 'skeletor'


if __name__ == '__main__':
    pytest.main()
