# Read an arbitrary h5 file in order to study its structure
# Goal: reveal and duplicate a CANSAS file such that it can be imported by SASVIEW
import numpy as np
from enum import Enum
import h5py


def parse_h5_entry(h5_entry):
    """Parse an HDF5 entry and generate an HDFNode object including all the sub entries

    Parameters
    ----------
    h5_entry: ~h5py._hl.dataset.Dataset, ~h5py._hl.group.Group, ~h5py._hl.files.File
        h5py entries including data set, group and file

    Returns
    -------
    HDFNode
        an HDFNode

    """
    # Create entry node instance
    entry_node = None
    for h5_entry_type, buffer_node_class in [(h5py._hl.files.File, FileNode),
                                             (h5py._hl.group.Group, GroupNode),
                                             (h5py._hl.files.File, DataSetNode)]:
        if isinstance(h5_entry, h5_entry_type):
            # generate node
            entry_node = buffer_node_class()
            # parse
            entry_node.parse_h5_entry(h5_entry)
            break
    # Check
    if entry_node is None:
        raise RuntimeError('HDF entry of type {} is not supported'.format(type(h5_entry)))

    return entry_node


class HDFNode(object):
    """
    an HDF node with more information
    """
    def __init__(self, name=None):
        """initialization

        Parameters
        ----------
        name: str, None
            entry name
        """
        self._name = name

        # Set attributes and etc
        self._attributes = dict()

        return

    def parse_h5_entry(self, h5_entry):
        """Parse an HDF5 entry

        Parameters
        ----------
        h5_entry

        Returns
        -------

        """
        # Parse attributes
        # Reset data structure
        self._attributes = dict()

        # parse attributes
        for attr_name in h5_entry.attrs:
            # returned h5_attribute in fact is attribute name
            self._attributes[attr_name] = h5_entry.attrs[attr_name]

        return

    @property
    def name(self):
        return self._name

    def write(self, input):
        """

        Parameters
        ----------
        input: str, ~h5py._hl.group.Group, ~h5py._hl.files.File
            Node to input

        Returns
        -------

        """
        raise NotImplementedError('Virtual method to write {}'.format(input))

    def write_attributes(self, curr_entry):

        # attributes
        for attr_name in self._attributes:
            curr_entry.attrs[attr_name] = self._attributes[attr_name]


class GroupNode(HDFNode):
    """
    Node for an HDF Group
    """
    def __init__(self):
        """
        Initialization
        """
        super(GroupNode, self).__init__()

        self._children = list()

    def parse_h5_entry(self, h5_entry):
        """Parse HDF5 entry

        Parameters
        ----------
        h5_entry: ~h5py._hl.group.Group
            hdf5 entry

        Returns
        -------
        None

        """
        # Parse in general way
        super(GroupNode, self).parse_h5_entry(h5_entry)

        # parse children
        children_names = h5_entry.keys()
        for child_name in children_names:
            child_node = parse_h5_entry(h5_entry[child_name])
            self._children.append(child_node)

    def write(self, parent_entry):
        """Write buffer node to an HDF entry

        Parameters
        ----------
        parent_entry:  ~h5py._hl.dataset.Dataset, ~h5py._hl.group.Group, ~h5py._hl.files.File
            parent HDF node

        Returns
        -------

        """
        # create group or data set
        # h5py._hl.group.Group only
        curr_entry = parent_entry.create_group(self._name)
        # write
        self.write_attributes(curr_entry)

    def write_content(self, curr_entry):
        # write child
        for child in self._children:
            child.write(curr_entry)

        # attributes
        self.write_attributes(curr_entry)


class FileNode(GroupNode):
    """
    Node for an HDF file
    """
    def __init__(self):
        """
        Initialization
        """
        super(FileNode, self).__init__()

    def write(self, file_name):
        """Write to a file

        Parameters
        ----------
        file_name: str
            Name of file to write to

        Returns
        -------

        """
        # create file node
        curr_entry = h5py.File(file_name, 'w')
        # write
        self.write_content(curr_entry)
        # close
        curr_entry.close()


class DataSetNode(HDFNode):
    """
    Node for data set
    """
    def __init__(self):
        """
        Initialization
        """
        super(DataSetNode, self).__init__()

        self._value = None

    def parse_h5_entry(self, h5_entry):
        """Parse HDF5 entry

        Parameters
        ----------
        h5_entry: ~h5py._hl.group.Group
            hdf5 entry

        Returns
        -------
        None

        """
        # Parse in general way
        super(DataSetNode, self).parse_h5_entry(h5_entry)

        # Parse value
        self._value = h5_entry.value

    def write(self, parent_entry):
        """Write buffer node to an HDF entry

        Parameters
        ----------
        parent_entry:  ~h5py._hl.group.Group, ~h5py._hl.files.File
            parent HDF node

        Returns
        -------

        """
        # Generate current entry and set the data
        curr_entry = parent_entry.create_dataset(self._name, data=self._value)

        self.write_attributes(curr_entry)
