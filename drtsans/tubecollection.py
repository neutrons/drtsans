import functools
from inspect import signature
import numpy as np
import re

from mantid.kernel import V3D
from mantid.api import mtd


class ElementComponentInfo:
    def __init__(self, component_info, component_info_index):
        r"""Wrapper to one of the components in ~mantid.ExperimentInfo.ComponentInfo

        Parameters
        ----------
        component_info: ~mantid.ExperimentInfo.ComponentInfo
        component_info_index: int
            Index to be used with the component_info object
        """
        self._component_info = component_info
        self.component_info_index = component_info_index

    def _decrement_arity(self, attribute, alternate_index=None):
        index = self.component_info_index if alternate_index is None else alternate_index
        if callable(attribute) is True:
            try:
                parameters = signature(attribute).parameters
            except ValueError:  # no signature found for builtin <Boost.Python.function object at ...>
                # Crappy hack. Determine number of arguments from docstring ! :(
                doc_string = attribute.__doc__
                n_parenthesis, first_index, last_index = 1, 1 + doc_string.find('('), len(doc_string)
                counter_incremental = {'(': 1, ')': -1}
                for i in range(first_index, len(doc_string)):
                    n_parenthesis += counter_incremental.get(doc_string[i], 0)
                    if n_parenthesis == 0:
                        last_index = i
                        break
                parameters = re.findall(',', doc_string[first_index: last_index])[0: -1]
            if len(parameters) == 0:
                return attribute(index)
            return functools.partial(attribute, index)
        return attribute

    def __getattr__(self, item):
        r"""
        Overload methods of ~mantid.ExperimentInfo.ComponentInfo so that they act as methods of ElementComponentInfo.

        Methods of ~mantid.ExperimentInfo.ComponentInfo with only one argument, the component info index, become
        read-only properties of ElementComponentInfo. Example: ``component_info.children(90177)``
         becomes ``element_component_info.children``

        Methods of ~mantid.ExperimentInfo.ComponentInfo with more than one argument where the first argument is the
        component info index become methods of ElementComponentInfo with this argument removed. Example:
        ``component_info.setPosition(90177, V3D(0.0, 0.0, 0.0)) becomes
        ``element_component_info.setPosition(V3D(0.0, 0.0, 0.0))``
        """
        _component_info = self.__dict__['_component_info']
        try:
            attribute = getattr(_component_info, item)
            return self._decrement_arity(attribute)
        except AttributeError:
            return getattr(self, item)


class PixelInfo(ElementComponentInfo):
    def __init__(self, component_info, component_info_index, detector_info=None):
        r"""Wrapper of ~mantid.ExperimentInfo.ComponentInfo when the component is a detector pixel."""
        super().__init__(component_info, component_info_index)
        self._detector_info = detector_info
        # Attributes pertaining to spectrumInfo. Maybe better if separate on another class
        self._spectrum_info = None
        self.spectrum_index = None

    def __getattr__(self, item):
        r"""
        Overload methods of ~mantid.ExperimentInfo.ComponentInfo and ~mantid.ExperimentInfo.DetectorInfo so that they
        act as methods of ElementComponentInfo.

        Methods of ~mantid.ExperimentInfo.ComponentInfo with only one argument, the component info index, become
        read-only properties of ElementComponentInfo. Example: ``component_info.children(90177)``
         becomes ``element_component_info.children``

        Methods of ~mantid.ExperimentInfo.ComponentInfo with more than one argument where the first argument is the
        component info index become methods of ElementComponentInfo with this argument removed. Example:
        ``component_info.setPosition(90177, V3D(0.0, 0.0, 0.0)) becomes
        ``element_component_info.setPosition(V3D(0.0, 0.0, 0.0))``
        """
        try:  # try method of componentInfo
            _component_info = self.__dict__['_component_info']
            attribute = getattr(_component_info, item)
            return self._decrement_arity(attribute)
        except AttributeError:
            try:  # try method of detectorInfo
                _detector_info = self.__dict__['_detector_info']
                attribute = getattr(_detector_info, item)
                return self._decrement_arity(attribute)
            except AttributeError:
                try:  # try method of spectrumInfo
                    _spectrum_info = self.__dict__['_spectrum_info']
                    attribute = getattr(_spectrum_info, item)
                    return self._decrement_arity(attribute, alternate_index=self.__dict__['spectrum_index'])
                except AttributeError:
                    return getattr(self, item)

    @property
    def detector_info(self):
        return self._detector_info

    @detector_info.setter
    def detector_info(self, det_info):
        if self._detector_info is not None:
            raise AttributeError('Detector info can be initialized only once')
        self._detector_info = det_info

    @property
    def position(self):
        return np.array(self._component_info.position(self.component_info_index))

    @position.setter
    def position(self, xyz):
        r"""Update the coordinates of the pixel

        Parameters
        ----------
        xyz: list
            Either a:
            - three-item iterable with the new X, Y, and Z coordinates
            - two-item tuple of the form ('x', float), or ('y', float), or ('z', float) if we only want to update one
              coordinate.
        """
        new_position = self.position
        if len(xyz) == 3:
            new_position = xyz  # substitute with new position
        else:
            new_position['xyz'.find(xyz[0])] = xyz[1]  # update selected coordinate
        self.setPosition(V3D(*list(new_position)))

    @property
    def width(self):
        return self.scaleFactor[0] * self.shape.getBoundingBox().width()[0]

    @width.setter
    def width(self, w):
        scale_factor = self.scaleFactor
        scale_factor[0] = w / self.shape.getBoundingBox().width()[0]
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def height(self):
        return self.scaleFactor[1] * self.shape.getBoundingBox().width()[1]

    @height.setter
    def height(self, h):
        scale_factor = self.scaleFactor
        scale_factor[1] = h / self.shape.getBoundingBox().width()[1]
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def area(self):
        return self.width * self.height

    #######
    #  Methods pertaining to spectrumInfo. Maybe better if separate on another class
    #######
    def insert_spectrum(self, spectrum_info, spectrum_index):
        if self._spectrum_info is not None:
            raise ValueError('Spectrum Info has already been inserted')
        self._spectrum_info = spectrum_info
        self.spectrum_index = spectrum_index


class TubeInfo(ElementComponentInfo):

    @staticmethod
    def is_valid_tube(component_info, component_index):
        children_indexes = [int(i) for i in component_info.children(component_index)]
        if len(children_indexes) == 0:
            return False  # we hit a detector
        if component_info.isDetector(children_indexes[0]) is False:
            return False  # quick check. The first child is not a detector
        children_are_detectors = [component_info.isDetector(index) for index in children_indexes]
        if len(set(children_are_detectors)) != 1:
            return False  # at least one child is not a detector
        return True

    def __init__(self, component_info, component_info_index):
        r"""Wrapper of ~mantid.ExperimentInfo.ComponentInfo when the component is a tube of detector pixels."""
        super().__init__(component_info, component_info_index)
        self._pixels = None
        if self.is_valid_tube(component_info, component_info_index) is False:
            raise ValueError('The component index is not associated to a valid tube')

    @property
    def pixels(self):
        if self._pixels is None:
            self._pixels = [PixelInfo(self._component_info, int(i)) for i in self.children]
        return self._pixels

    def __getitem__(self, item):
        return self.pixels[item]  # iterate over the pixels


class TubeCollection(ElementComponentInfo):

    @staticmethod
    def map_detector_to_spectrum(input_workspace):
        r"""A map from detectorInfo index (or componentInfo index) to workspace spectrum index"""
        spectrum_info = input_workspace.spectrumInfo()
        detector_to_spectrum = dict()
        for spectrum_index in range(input_workspace.getNumberHistograms()):
            detector_info_index = spectrum_info.getSpectrumDefinition(spectrum_index)[0][0]
            detector_to_spectrum[detector_info_index] = spectrum_index
        return detector_to_spectrum

    def __init__(self, input_workspace, component_name):
        workspace_handle = mtd[str(input_workspace)]
        component_info = workspace_handle.componentInfo()
        for component_index in range(component_info.root(), -1, -1):
            if component_info.name(component_index) == component_name:
                super().__init__(component_info, component_index)
                break
            if component_info.isDetector(component_index) is True:  # we reached the detector with the highest index
                raise RuntimeError(f'Could not find a component with name "{component_name}"')
        self._tubes = list()
        self._sorting_permutations = {}
        self._input_workspace = workspace_handle
        # A map from detectorInfo index (or componentInfo index) to workspace spectrum index
        self.detector_to_spectrum = self.map_detector_to_spectrum(workspace_handle)

    def __getitem__(self, item):
        return self.tubes[item]

    @property
    def tubes(self):
        r"""List of TubeInfo objects ordered using their component index, from smallest to highest index."""
        if len(self._tubes) == 0:
            # Find the mapping between spectrum indexes and detectorInfo indexes
            spectrum_info = self._input_workspace.spectrumInfo()
            detector_info = self._input_workspace.detectorInfo()
            non_detector_indexes = sorted([int(i) for i in set(self.componentsInSubtree)-set(self.detectorsInSubtree)])
            for component_index in non_detector_indexes:
                try:
                    tube = TubeInfo(self._component_info, component_index)
                except ValueError:  # the component index is not associated to a tube
                    continue
                for pixel in tube:
                    pixel.detector_info = detector_info
                    pixel.insert_spectrum(spectrum_info, self.detector_to_spectrum[pixel.component_info_index])
                self._tubes.append(tube)
        return self._tubes

    @property
    def sorted_views(self):
        return self._sorting_permutations.keys()

    def sorted(self, key=None, reverse=False, view='decreasing X'):
        r"""
        List of tubes in the prescribed order

        Parameters
        ----------
        key: :py:obj:`function`
            Function of one argument to extract a comparison key from each element in ``self.tubes``. If None, then
            the selected ``view`` determines the order
        reverse: bool
            Reverse the order resulting from application of ``key`` or ``view``
        view: str
            Built-in permutations of the tubes prescribing a particular order. Valid views are:
            - 'decreasing X' order the tubes by decreasing X-coordinate. This view flattens a double detector panel and
              is viewed 'from left to right' when viewed from the sample.
            - 'spectrum index' order the tubes by increasing spectrum index (workspace index)
        """
        if key is not None:
            return sorted(self._tubes, key=key, reverse=reverse)
        permutation = self._sorting_permutations.get(view, None)
        if permutation is None:
            if view == 'decreasing X':  # initialize this view
                x_coords = [tube[0].position[0] for tube in self.tubes]  # X coords for first pixel of each tube
                permutation = np.flip(np.argsort(x_coords), axis=0).tolist()
                self._sorting_permutations['decreasing X'] = permutation
            elif view == 'spectrum index':  # initialize this view
                # spectrum index of first pixel for each tube
                permutation = np.argsort([tube[0].spectrum_index for tube in self.tubes])
                self._sorting_permutations['spectrum index'] = permutation
        sorted_list = [self._tubes[i] for i in permutation]
        return sorted_list if reverse is False else sorted_list[::-1]
