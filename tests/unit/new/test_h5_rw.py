import pytest
from drtsans.h5_buffer import HDFNode, GroupNode, DataSetNode


def test_create_base_node():
    """

    Returns
    -------

    """
    attributes = {'A': 3, 'B': 'hello world'}
    attributes_alt1 = {'A': 3, 'B': 'hello worlds'}
    attributes_alt2 = {'A': 3, 'B': 'hello world', 'C': 3.23}
    name = 'test'

    node1 = HDFNode(name=name)
    node1.add_attributes(attributes)

    node2 = HDFNode(name=name)
    node2.add_attributes(attributes)

    node3 = HDFNode(name=name)
    node3.add_attributes(attributes_alt1)

    node4 = HDFNode(name=name)
    node4.add_attributes(attributes_alt2)

    # node 1 and 2 shall be same
    node1.match(node2)

    # node 1 and 3 shall have Value error
    with pytest.raises(ValueError):  # , 'Expecting a ValueError'):
        node1.match(node3)

    # node 1 and 4 shall have KeyError
    with pytest.raises(KeyError):
        node1.match(node4)


if __name__ == '__main__':
    pytest.main(__file__)
