import functools
from inspect import signature
import numpy as np

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

    def _decorate_component_method(self, attribute):
        if callable(attribute) is True:
            if len(signature(attribute).parameters) == 0:
                return attribute(self._component_info_index)
            return functools.partial(attribute, self._component_info_index)
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
            return self._decorate_component_method(attribute)
        except AttributeError:
            return getattr(self, item)


class PixelInfo(ElementComponentInfo):
    def __init__(self, component_info, component_info_index, detector_info=None):
        r"""Wrapper of ~mantid.ExperimentInfo.ComponentInfo when the component is a detector pixel."""
        self._detector_info = detector_info
        super(self).__init__(component_info, component_info_index)

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
        try:
            _component_info = self.__dict__['_component_info']
            attribute = getattr(_component_info, item)
            return self._decorate_component_method(attribute)
        except AttributeError:
            try:
                _detector_info = self.__dict__['_detector_info']
                attribute = getattr(_detector_info, item)
                return self._decorate_component_method(attribute)
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
        return self._component_info(self._component_info_index).position()

    @position.setter
    def position(self, xyz=None, x=None, y=None, z=None):
        r"""Update the coordinates of the pixel

        Parameters
        ----------
        xyz: list
            A 3-number iterable with the new coordinates
        x, y, z: float
            Update only coordinates along these axes"""
        new_position = self.position
        if xyz is not None:
            new_position = xyz
        for coord, index in zip((x, y, z), (0, 1, 2)):
            if coord is not None:
                new_position[index] = coord
        self.setPosition(V3D(*list(new_position)))

    @property
    def width(self):
        return self.scaleFactor(0) * self.getBoundingBox.width(0)

    @width.setter
    def width(self, w):
        scale_factor = self.scaleFactor
        scale_factor[0] = w / self.getBoundingBox.width(0)
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def height(self):
        return self.scaleFactor(1) * self.getBoundingBox.width(1)

    @height.setter
    def height(self, h):
        scale_factor = self.scaleFactor
        scale_factor[1] = h / self.getBoundingBox.width(1)
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def area(self):
        return self.width * self.height


class TubeInfo(ElementComponentInfo):

    @staticmethod
    def is_valid_tube(cls, component_info, component_index):
        children_indexes = component_info.children(component_index)
        if component_info.isDetector(children_indexes[0]) is False:
            return False  # quick check. The first child is not a detector
        children_detector = [component_info.isDetector(index) for index in children_indexes]
        if len(set(children_detector)) != 1:
            return False  # at least one child is not a detector
        return True

    def __init__(self, component_info, component_info_index):
        r"""Wrapper of ~mantid.ExperimentInfo.ComponentInfo when the component is a tube of detector pixels."""
        self._pixels = None
        if self.is_valid_tube(component_info, component_info_index) is False:
            raise ValueError('The component index is not associated to a valid tube')
        super(self).__init__(component_info, component_info_index)

    @property
    def pixels(self):
        if self._pixels is None:
            component_indexes = self._component_info.children(self._component_info_index)
            self._pixels = [PixelInfo(self._component_info, component_index) for component_index in component_indexes]
        return self._pixels

    def __getitem__(self, item):
        return self.pixels[item]  # iterate over the pixels


class TubeCollection(ElementComponentInfo):

    def __init__(self, input_workspace, component_name):
        self._tubes = None
        self._sorting_permutations = {}
        self.input_workspace = mtd[str(input_workspace)]
        component_info = self._input_workspace.componentInfo()
        for component_index in range(self.component_info.root(), -1, -1):
            if component_info.name(component_index) == component_name:
                super(self).__init__(component_info, component_index)
                break
            if component_info.isDetector(component_index) is True:
                raise RuntimeError(f'Could not find a component with name "{component_name}"')

    def __getitem__(self, item):
        return self.tubes[item]

    @property
    def tubes(self):
        r"""List of TubeInfo objects ordered using their component index, from smallest to highest index."""
        if self._tubes is None:
            detector_info = self._input_workspace.detectorInfo()
            for component_index in sorted(list(set(self.componentsInSubtree) - set(self.detectorsInSubtree))):
                try:
                    tube = TubeInfo(self._component_info, component_index)
                    for pixel in tube:
                        pixel.detector_info = detector_info
                    self._tubes.append(tube)
                except ValueError:  # the component index is not associated to a tube
                    continue
        return self._tubes

    @property
    def sorted_views(self):
        return self._sorting_permutations.keys()

    def sorted(self, key=None, reverse=False, view='increasing X'):
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
            'increasing X' order the tubes by increasing X-coordinate. This view flattens a double detector panel.
        """
        if key is not None:
            return sorted(self._tubes, key=key, reverse=reverse)
        permutation = self._sorting_permutations.get(view, None)
        if permutation is None:
            if view == 'increasing X':  # initialize this view
                x_coords = [tube.position[0] for tube in self.tubes]
                self._sorting_permutations['increasing X'] = np.argsort(x_coords).tolist()
        sorted_list = [self._tubes[i] for i in permutation]
        return sorted_list if reverse is False else sorted_list[::-1]
