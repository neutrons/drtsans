import functools
from inspect import signature
import numpy as np

from mantid.kernel import V3D
from mantid.api import mtd


def _decrement_arity(attribute, index):
    r"""
    Decrement the arity of callable ``attribute`` by either evaluating this function, or creating a partial
    function.

    If the arity of ``attribute`` is one, then we evaluate ``attribute`` by passing ``index`` as its argument,
    finally returning the value resulting from the evaluation. If the arity is bigger than one, then we create a
    partial function by passing ``index`` as the first argument of ``attribute``.

    Parameters
    ----------
    attribute: function
        Callable to evaluate or reduce it's arity with a partial function.
    index: int
        Value to use as fist argument of callable ``attribute``.

    Returns
    -------
    PyObject, function
    """
    if callable(attribute) is True:
        try:
            parameters = list(signature(attribute).parameters)
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
            parameters = doc_string[first_index: last_index].split(',')
            if parameters[0][-4:] == 'self':  # `attribute` is a bound method
                parameters.pop()
        if len(parameters) == 0:
            result = attribute()
            return np.array(result) if isinstance(result, V3D) else result
        if len(parameters) == 1:
            result = attribute(index)
            return np.array(result) if isinstance(result, V3D) else result
        return functools.partial(attribute, index)
    return attribute  # nothing to do is `attribute` is not callable


def _inverse_map(list_of_functions, *args, **kwargs):
    r"""Apply a list of functions to a set of given arguments"""
    return [function(*args, **kwargs) for function in list_of_functions]


class SpectrumInfo:

    def __init__(self, input_workspace, workspace_index):
        r"""Wrapper to ~mantid.api.SpectrumInfo. We reduce the arity for the methods of `SpectrumInfo`, thus
        converting them into methods of `SpectrumInfo`.

        The class contains additional methods for SpectrumInfo that are wrappers to methods
        of ~mantid.api.Workspace such as `readY`

        Parameters
        ----------
        input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        workspace_index: int, list
            Single index or list of indexes to be used with the ``spectrum_info`` object.
        """
        self._workspace = mtd[str(input_workspace)]
        self._spectrum_info = self._workspace.spectrumInfo()
        self.spectrum_info_index = workspace_index

    def __getattr__(self, item):
        r"""
        Overload methods of ~mantid.api.SpectrumInfo so that they act as methods of SpectrumInfo.

        Methods of ~mantid.api.SpectrumInfo with only one argument, the component info index, become
        read-only properties of SpectrumInfo. Example: ``spectrum_info.hasUniqueDetectors(90177)``
         becomes ``element_component_info.hasUniqueDetectors``

        Methods of ~mantid.geometry.SpectrumInfo with more than one argument where the first argument is the
        component info index become methods of SpectrumInfo with this argument removed.

        Parameters
        ----------
        item: str
        """
        _spectrum_info = self.__dict__['_spectrum_info']
        try:
            attribute = getattr(_spectrum_info, item)
            if isinstance(self.spectrum_info_index, int):
                return _decrement_arity(attribute, self.spectrum_info_index)
            # `self.spectrum_info_index` is a list of spectrum indexes
            arity_decremented_attributes = [_decrement_arity(attribute, i) for i in self.spectrum_info_index]
            if callable(arity_decremented_attributes[0]) is True:  # attribute's arity was bigger than one
                # A partial function that will call each arity-decremented attribute with remaining method arguments
                return functools.partial(_inverse_map, arity_decremented_attributes)
            return arity_decremented_attributes  # attribute's arity was less than two
        except AttributeError:
            return super().__getattr__(item)  # Next class in the Method Resolution Order

    def _iterate_with_indexes(self, function_name, array_type=np.array):
        r"""
        Utility function to discriminate between one or more indexes contained in ``spectrum_info_index`` when
        evaluating a method invoked via a Mantid's workspace object.

        Parameters
        ----------
        function_name: str
            Method name invoked with the workspace handle ``self._workspace``
        array_type: type
            Data structure constructor when ``self.spectrum_info_index`` is a list of indexes. Examples: list,
            numpy.array

        Returns
        -------
        Object
            Object returned by invoking the function over the indexes in ``self.spectrum_info_index``.
        """
        function = getattr(self._workspace, function_name)
        if isinstance(self.spectrum_info_index, int):
            return function(self.spectrum_info_index)
        return array_type([function(index) for index in self.spectrum_info_index])

    @property
    def readX(self):
        return self._iterate_with_indexes('readX')

    @property
    def readY(self):
        return self._iterate_with_indexes('readY')

    @property
    def readE(self):
        return self._iterate_with_indexes('readE')

    def __len__(self):
        return len(self.children)


class ElementComponentInfo:
    def __init__(self, component_info, component_info_index):
        r"""Wrapper to one of the components in ~mantid.ExperimentInfo.ComponentInfo

        Parameters
        ----------
        component_info: ~mantid.geometry.ComponentInfo
        component_info_index: int
            Index to be used with the component_info object
        """
        self._component_info = component_info
        self.component_info_index = component_info_index

    def __getattr__(self, item):
        r"""
        Overload methods of ~mantid.geometry.ComponentInfo so that they act as methods of ElementComponentInfo.

        Methods of ~mantid.geometry.ComponentInfo with only one argument, the component info index, become
        read-only properties of ElementComponentInfo. Example: ``component_info.children(90177)``
         becomes ``element_component_info.children``

        Methods of ~mantid.geometry.ComponentInfo with more than one argument where the first argument is the
        component info index become methods of ElementComponentInfo with this argument removed. Example:
        ``component_info.setPosition(90177, V3D(0.0, 0.0, 0.0)) becomes
        ``element_component_info.setPosition(V3D(0.0, 0.0, 0.0))``
        """
        _component_info = self.__dict__['_component_info']
        try:
            attribute = getattr(_component_info, item)
            return _decrement_arity(attribute, self.component_info_index)
        except AttributeError:
            return super().__getattr__(item)  # Next class in the Method Resolution Order

    @property
    def children(self):
        return [int(index) for index in self._component_info.children(self.component_info_index)]  # cast to int


def _resolve_indexes(input_workspace, component_info_index, workspace_index):
    r"""Resolve the component info index or the workspace index when only one of the two indexes is provided"""
    if workspace_index is not None and component_info_index is not None:
        return component_info_index, workspace_index  # nothing to do if both were provided
    if workspace_index is None and component_info_index is None:
        raise RuntimeError('Either "component/detector_info_index" or "spectrum_info/workspace_index" must be passed')

    input_workspace = mtd[str(input_workspace)]
    get_spectrum_definition = input_workspace.spectrumInfo().getSpectrumDefinition

    if component_info_index is None:
        component_info_index = get_spectrum_definition(workspace_index)[0][0]

    if workspace_index is None:
        def get_component_index(index):
            return get_spectrum_definition(index)[0][0]
        # Expensive search, since no a priori sorting between workspace indexes and component info indexes!!
        for index in range(0, input_workspace.getNumberHistograms()):
            if get_component_index(index) == component_info_index:
                workspace_index = index
                break

    return component_info_index, workspace_index


class PixelSpectrum(ElementComponentInfo, SpectrumInfo):
    def __init__(self, input_workspace, component_info_index=None, workspace_index=None):
        r"""
        Wrapper of ~mantid.geometry.ComponentInfo, ~mantid.api.DetectorInfo, and ~mantid.api.SpectrumInfo for a pixel.

        Additionally, the class contains additional methods invoked through objects of type ~mantid.api.Workspace
        such as `readY`.

        Parameters
        ----------
        input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        component_info_index: int
            Index to be used with the ``ElementComponentInfo`` object. If not supplied, then ``workspace_index``
            is used for initialization.
        workspace_index: int
            If not supplied, then ``detector_info_index`` is used for initialization.
        """
        component_info_index, workspace_index = _resolve_indexes(input_workspace, component_info_index,
                                                                 workspace_index)
        input_workspace = mtd[str(input_workspace)]
        self._detector_info = input_workspace.detectorInfo()
        ElementComponentInfo.__init__(self, input_workspace.componentInfo(), component_info_index)
        SpectrumInfo.__init__(self, input_workspace, workspace_index)

    def __getattr__(self, item):
        try:
            _detector_info = self.__dict__['_detector_info']
            attribute = getattr(_detector_info, item)
            return _decrement_arity(attribute, self.component_info_index)
        except AttributeError:
            return super().__getattr__(item)  # next class in the Method Resolution Order, (SpectrumInfo.__getattr__)

    @property
    def position(self):
        r"""Cartesian coordinates of the pixel detector.

        Coordinates typically correspond to the center of a cuboid detector, or the center for the base of a
        cylindrical pixel.
        """
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
              of the three coordinates.
        """
        new_position = self.position
        if len(xyz) == 3:
            new_position = xyz  # substitute with new position
        else:
            new_position['xyz'.find(xyz[0])] = xyz[1]  # update selected coordinate
        self.setPosition(V3D(*list(new_position)))

    @property
    def width(self):
        r"""Extent of the pixel detector along the Y-axis"""
        return self.scaleFactor[0] * self.shape.getBoundingBox().width()[0]

    @width.setter
    def width(self, w):
        scale_factor = self.scaleFactor
        scale_factor[0] = w / self.shape.getBoundingBox().width()[0]
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def height(self):
        r"""Extent of the pixel detector along the Z-axis"""
        return self.scaleFactor[1] * self.shape.getBoundingBox().width()[1]

    @height.setter
    def height(self, h):
        scale_factor = self.scaleFactor
        scale_factor[1] = h / self.shape.getBoundingBox().width()[1]
        self.setScaleFactor(V3D(*scale_factor))

    @property
    def area(self):
        r"""Product of pixel width and height"""
        return self.width * self.height


class TubeSpectrum(ElementComponentInfo, SpectrumInfo):

    @staticmethod
    def is_valid_tube(component_info, component_index):
        r"""
        Determine if all the immediate children-components are detectors, as in the case of a tube.

        Warnings
        --------
        This function does not invalidate other components than tubes whose immediate children-components are
        all detectors.

        Parameters
        ----------
        component_info: ~mantid.geometry.ComponentInfo
        component_index: int

        Returns
        -------
        bool
        """
        children_indexes = [int(i) for i in component_info.children(component_index)]
        if len(children_indexes) == 0:
            return False  # we hit a detector
        if component_info.isDetector(children_indexes[0]) is False:
            return False  # quick check. The first child is not a detector
        children_are_detectors = [component_info.isDetector(index) for index in children_indexes]
        if len(set(children_are_detectors)) != 1:
            return False  # at least one child is not a detector
        return True

    def __init__(self, input_workspace, component_info_index, workspace_indexes):
        r"""Wrapper of ~mantid.geometry.ComponentInfo when the component is a tube of detector pixels.

        Parameters
        ----------
        input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        component_info_index: int
            Index corresponding to the tube
        """
        input_workspace = mtd[str(input_workspace)]
        if self.is_valid_tube(input_workspace.componentInfo(), component_info_index) is False:
            raise ValueError('The component index is not associated to a valid tube')
        self._pixels = list()
        SpectrumInfo.__init__(self, input_workspace, workspace_indexes)
        ElementComponentInfo.__init__(self, input_workspace.componentInfo(), component_info_index)

    @property
    def pixels(self):
        r"""List of ~drtsans.tubecollection.PixelSpectrum objects making up the tube.

        Returns
        -------
        list
        """
        if len(self._pixels) == 0:
            for component_info_index, workspace_index in zip(self.children, self.spectrum_info_index):
                self._pixels.append(PixelSpectrum(self._workspace, component_info_index=component_info_index,
                                                  workspace_index=workspace_index))
        return self._pixels

    def __getitem__(self, item):
        return self.pixels[item]  # iterate over the pixels

    def __len__(self):
        return len(self.spectrum_info_index)


class TubeCollection(ElementComponentInfo):

    @staticmethod
    def map_detector_to_spectrum(input_workspace):
        r"""A mapping from detector info index (or component info index) to spectrum index (or workspace index)

        Parameters
        ----------
        input_workspace: str, ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace

        Returns
        -------
        dict
        """
        input_workspace = mtd[str(input_workspace)]
        get_spectrum_definition = input_workspace.spectrumInfo().getSpectrumDefinition

        def get_detector_info_index(workspace_index):
            return get_spectrum_definition(workspace_index)[0][0]

        detector_to_spectrum = dict()
        for spectrum_index in range(input_workspace.getNumberHistograms()):
            detector_to_spectrum[get_detector_info_index(spectrum_index)] = spectrum_index

        return detector_to_spectrum

    def __init__(self, input_workspace, component_name):
        r"""
        A list of ~drtsans.tubecollection.TubeInfo objects making up one of the instrument components, such as the
        front panel, a double panel, or the whole instrument.

        Parameters
        ----------
        input_workspace: ~mantid.api.MatrixWorkspace, ~mantid.api.IEventWorkspace
        component_name: str
            One of the named components in the instrument geometry file.
        """
        input_workspace = mtd[str(input_workspace)]
        component_info = input_workspace.componentInfo()
        for component_index in range(component_info.root(), -1, -1):
            if component_info.name(component_index) == component_name:
                super().__init__(component_info, component_index)
                break
            if component_info.isDetector(component_index) is True:  # we reached the detector with the highest index
                raise RuntimeError(f'Could not find a component with name "{component_name}"')
        self._tubes = list()
        self._sorting_permutations = {}
        self._input_workspace = input_workspace
        # A map from detectorInfo index (or componentInfo index) to workspace spectrum index
        self.detector_to_spectrum = None

    def __getitem__(self, item):
        return self.tubes[item]

    def __len__(self):
        return len(self.tubes)

    @property
    def tubes(self):
        r"""List of ~drtsans.tubecollection.TubeSpectrum objects ordered using their component info indexes,
        from smallest to highest index."""
        if len(self._tubes) == 0:
            # Initialize the mapping between spectrum indexes and detectorInfo indexes
            self.detector_to_spectrum = self.map_detector_to_spectrum(self._input_workspace)
            # Iterate over the components of the instrument that are not detectors
            non_detector_indexes = sorted([int(i) for i in set(self.componentsInSubtree)-set(self.detectorsInSubtree)])
            for component_info_index in non_detector_indexes:
                if TubeSpectrum.is_valid_tube(self._component_info, component_info_index) is True:
                    tube_info = ElementComponentInfo(self._component_info, component_info_index)
                    # Find workspace indexes associated to the component/detector info indexes
                    workspace_indexes = [self.detector_to_spectrum[index] for index in tube_info.children]
                    self._tubes.append(TubeSpectrum(self._input_workspace, component_info_index, workspace_indexes))
        return self._tubes

    @property
    def sorted_views(self):
        r"""Dictionary mapping a particular sorting of the tubes to the permutation of their component info indexes
        required to attaing the desired sorting."""
        return self._sorting_permutations.keys()

    def sorted(self, key=None, reverse=False, view='decreasing X'):
        r"""
        Sort the list of ~drtsans.tubecollection.TubeInfo objects in the prescribed order.

        Parameters
        ----------
        key: :py:obj:`function`
            Function of one argument to extract a comparison key from each element in ``self.tubes``. If None, then
            the selected ``view`` determines the order
        reverse: bool
            Reverse the order resulting from application of ``key`` or ``view``
        view: str
            Built-in permutations of the tubes prescribing a particular order. Valid views are:
            - 'decreasing X': order the tubes by decreasing X-coordinate. This view can "flatten" a double
            detector panel when viewed from the sample "from left to right".
            - 'workspace index': order the tubes by increasing workspace index for the first pixel of each tube.
        """
        if key is not None:
            return sorted(self._tubes, key=key, reverse=reverse)
        permutation = self._sorting_permutations.get(view, None)
        if permutation is None:
            if view == 'decreasing X':  # initialize this view
                x_coords = [tube.position[0] for tube in self.tubes]  # X coords of each tube
                permutation = np.flip(np.argsort(x_coords), axis=0).tolist()
                self._sorting_permutations['decreasing X'] = permutation
            elif view == 'workspace index':  # initialize this view
                # spectrum index of first pixel for each tube
                permutation = np.argsort([tube[0].spectrum_info_index for tube in self.tubes])
                self._sorting_permutations['spectrum index'] = permutation
        sorted_list = [self._tubes[i] for i in permutation]
        return sorted_list if reverse is False else sorted_list[::-1]
