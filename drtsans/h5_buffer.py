# Read an arbitrary h5 file in order to study its structure
# Goal: reveal and duplicate a CANSAS file such that it can be imported by SASVIEW
import numpy as np
from enum import Enum


class HDFDictKeys(Enum):
    attributes = 'Attributes'
    children = 'SubNodes'
    data_value = 'Value'
    data_type = 'Type'


class HDFNode(object):
    """
    an HDF node with more information
    """
    def __init__(self, h5entry, parent):
        """initialization

        Parameters
        ----------
        h5entry: HDF entry
        parent: HDFNode, None
            for parent to create linked list
        """
        # print('Init: h5entry = {}'.format(h5entry))
        self._name = h5entry.name

        self._attributes, self._value, child_h5_entries = self.parse_h5node(h5entry)

        self._parent = parent

        # parse children
        self._children = list()
        for child_entry in child_h5_entries:
            child_node = HDFNode(h5entry[child_entry], self)
            self._children.append(child_node)

    @staticmethod
    def parse_h5node(h5_entry):
        """

        Parameters
        ----------
        h5_entry

        Returns
        -------
        tuple
            attributes (list), value (list), child hdf5 entry (list)

        """
        # Init data structure
        attr_dict = dict()

        # parse attributes
        for attr_name in h5_entry.attrs:
            # returned h5_attribute in fact is attribute name
            attr_dict[attr_name] = h5_entry.attrs[attr_name]

        # child or value
        try:
            # h5py._hl.group.Group
            children_names = h5_entry.keys()
            value = None
        except AttributeError:
            # h5py._hl.dataset.Dataset
            value = h5_entry[()]
            children_names = list()

        return attr_dict, value, children_names

    @property
    def name(self):
        return self._name

    def nice(self):
        """Print out nicely

        Returns
        -------

        """
        print('Name: {}'.format(self._name))
        print('Attributes: {}'.format(self._attributes.keys()))
        for child in self._children:
            print('Child: {}'.format(child.name))

        for child in self._children:
            child.nice()

    def write(self, parent_entry):
        """

        Returns
        -------

        """
        # print('Create entry: {}'.format(self._name))
        if self._name.endswith('_spectrum_sample'):
            print('SKIP!!!!!  Node {}'.format(self._name))
            return

        # create group or data set
        if self._value is None:
            # h5py._hl.group.Group
            if self._name != '/':
                curr_entry = parent_entry.create_group(self._name)
            else:
                # root entry is special
                curr_entry = parent_entry

            # write child
            for child in self._children:
                child.write(curr_entry)
        else:
            # h5py._hl.dataset.Dataset
            curr_entry = parent_entry.create_dataset(self._name, data=self._value)

        # attributes
        for attr_name in self._attributes:
            curr_entry.attrs[attr_name] = self._attributes[attr_name]

    def convert_to_dict(self):
        """Convert this node and the children to dictionary

        Returns
        -------

        """
        # convert node to dict entry
        node_dict = {HDFDictKeys.attributes.value: dict(),
                     HDFDictKeys.children.value: dict(),
                     HDFDictKeys.data_value.value: None,
                     HDFDictKeys.data_type.value: None}

        # attribute
        for attr_name in self._attributes:
            node_dict[HDFDictKeys.attributes.value][attr_name] = dict()
            node_dict[HDFDictKeys.attributes.value][attr_name][HDFDictKeys.data_value.value] =\
                str(self._attributes[attr_name])
            node_dict[HDFDictKeys.attributes.value][attr_name][HDFDictKeys.data_type.value] =\
                str(type(self._attributes[attr_name]))

        # child or value
        if self._value is None:
            # child: recursive
            for index, child in enumerate(self._children):

                # skip DASlogs but leave a few as template
                if self._name == 'DASlogs' and index > 1:
                    # skip but the first 3
                    continue
                # only keep bank1_events
                if self._name.startswith('bank'):
                    try:
                        bank_number = int(self._name.split('bank')[1].split('_')[0])
                        if bank_number > 1:
                            continue
                    except TypeError:
                        pass

                node_dict[HDFDictKeys.children.value][child.name] = child.convert_to_dict()
        else:
            # value
            if isinstance(self._value, np.ndarray):
                node_dict[HDFDictKeys.data_value.value] = str(self._value.shape)
                node_dict[HDFDictKeys.data_type.value] = str(self._value.dtype)
            else:
                node_dict[HDFDictKeys.data_value.value] = str(self._value)
                node_dict[HDFDictKeys.data_type.value] = str(type(self._value))

        return node_dict
