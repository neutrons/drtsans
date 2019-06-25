from __future__ import absolute_import, division, print_function


class Component:
    '''
    Stores information about the component
    '''

    # class variables that will cache the component details
    dim_x = -1
    dim_y = -1
    dims = -1
    first_index = -1

    def __init__(self, workspace, component_name):
        self._workspace = workspace
        self._initialization(component_name)

    def _initialization(self, component_name):
        first_det_id = self._detector_details(component_name)
        self._detector_first_ws_index(first_det_id)

    def _detector_details(self, component_name):
        """Private function that reads the instrument and get component_name
        dimensions and first detector id
        """
        instrument = self._workspace.getInstrument()
        component = instrument.getComponentByName(component_name)
        self.dim_x = component.nelements()

        if component[0].nelements() == 1:
            # Handles EQSANS
            self.dim_y = component[0][0].nelements()
            first_det_id = component[0][0][0].getID()
        else:
            # Handles BioSANS/GPSANS
            self.dim_y = component[0].nelements()
            first_det_id = component[0][0].getID()
        self.dims = self.dim_x * self.dim_y
        return first_det_id

    def _detector_first_ws_index(self, first_det_id):
        ''' sets the first_index of this component
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
