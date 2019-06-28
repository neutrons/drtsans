from __future__ import absolute_import, division, print_function


class Component:
    '''
    Stores information about the component
    '''

    # class variables that will cache the component details
    dim_x = -1
    dim_y = -1
    dims = -1  # total number of pixels
    first_index = -1  # workspace index of the smallest pixel id

    def __init__(self, workspace, component_name):
        self._workspace = workspace
        self._initialization(component_name)

    def _initialization(self, component_name):
        first_det_id = self._detector_details(component_name)
        self._detector_first_ws_index(first_det_id)

    def _num_pixels_in_tube(self, info, component_index):
        '''Recursive function that determines how many pixels are in a single
        tube (y-dimension). This assumes that things without grand-children
        are tubes'''
        component_index = int(component_index)
        children = info.children(component_index)

        grandchildren = info.children(int(children[0]))
        if len(grandchildren) == 0:
            return children[0], len(children)
        else:
            return self._num_pixels_in_tube(info, children[0])

    def _detector_details(self, component_name):
        """Private function that reads the instrument and get component_name
        dimensions and first detector id
        """
        component_info = self._workspace.componentInfo()
        detector_info = self._workspace.detectorInfo()
        component_index = component_info.indexOfAny(component_name)

        total_pixels = len(component_info.detectorsInSubtree(component_index))
        tube_index, self.dim_y = \
            self._num_pixels_in_tube(component_info, component_index)
        self.dim_x = total_pixels // self.dim_y
        self.dims = total_pixels

        return detector_info.detectorIDs()[tube_index]

    def _detector_first_ws_index(self, first_det_id):
        '''sets the first_index of this component
        '''
        for ws_index in range(self._workspace.getNumberHistograms()):
            if self._workspace.getSpectrum(
                    ws_index).getDetectorIDs()[0] == first_det_id:
                self.first_index = ws_index
                break
                return
        else:
            raise ValueError('Iterared WS and did not find first det id = '
                             '{}'.format(first_det_id))
