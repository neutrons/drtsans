import pytest
import numpy as np
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


def test_create_data_node():
    """Test creating and setting up DataNode

    Returns
    -------

    """
    attributes = {'A': 3, 'B': 'hello world'}
    attributes_alt1 = {'A': 3, 'B': 'hello worlds'}
    name = 'test data'

    node1 = DataSetNode(name)
    node1.add_attributes(attributes)
    node1.set_value(np.ndarray([1, 2, 3]))

    node2 = DataSetNode(name)
    node2.add_attributes(attributes)
    node2.set_value(np.ndarray([1, 2, 3]))

    node3 = DataSetNode(name)
    node3.add_attributes(attributes_alt1)
    node3.set_value(np.ndarray([1, 2, 3]))

    node4 = DataSetNode(name)
    node4.add_attributes(attributes)
    node4.set_value(np.ndarray([1, 2, 3.2]))

    # Node 1 and 2 shall be same
    node2.match(node1)

    # Node 1 and Node 3 are different with attributes
    with pytest.raises(ValueError):
        node1.match(node3)

    # Node 1 and Node 4 are different by value
    with pytest.raises(ValueError):
        node1.match(node4)


def test_check_type():
    """Test to check node type

    Returns
    -------

    """
    attributes = {'A': 3, 'B': 'hello world'}

    node1 = HDFNode('base')
    node1.add_attributes(attributes)

    node2 = GroupNode('group')
    node2.add_attributes(attributes)

    node3 = DataSetNode('data')
    node3.add_attributes(attributes)

    # node 3 shall match to node 1
    node1.match(node3)

    # node 1 does not match to node 3
    with pytest.raises(TypeError):
        node3.match(node1)

    # node 2 and node 3 are different
    with pytest.raises(TypeError):
        node2.match(node3)
    with pytest.raises(TypeError):
        node3.match(node2)


if __name__ == '__main__':
    pytest.main(__file__)
