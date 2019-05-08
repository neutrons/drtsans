

class Detector(object):

    _n_monitors = None

    _components_dimensions = {}

    def __init__(self, workspace):
        self._workspace = workspace
        self._instrument = workspace.getInstrument()

    def _get_detector_dimensions(self, component_name):
        """Private function that reads the instrument and get component_name
        dimensions
        """

        component = self._instrument.getComponentByName(component_name)
        n_tubes = component.nelements()

        if component[0].nelements() == 1:
            # Handles EQSANS
            n_pixels_per_tube = component[0][0].nelements()
        else:
            # Handles BioSANS/GPSANS
            n_pixels_per_tube = component[0].nelements()
        return n_tubes, n_pixels_per_tube

    def get_detector_dimensions(self, component_name):
        """ Reads the cached dimensions
        Parameters
        ----------
        component_name : string
            Valid component name

        Returns
        -------
        tuple
            tuple of integers x,y
        """

        if component_name in self._components_dimensions:
            return self._components_dimensions[component_name]
        else:
            x, y = self._get_detector_dimensions(component_name)
            self._components_dimensions[component_name] = (x, y)
            return x, y

    def get_number_of_monitors(self):
        """Get the number of monitors for this workspace.
        It is assumed that the monitors are in the beggining of the
        spectra list.

        Returns
        -------
        int
            The number of monitors
        """
        if not self._n_monitors:
            spectrum_info = self._workspace.spectrumInfo()
            monitor_index_list = []
            for ws_index in range(self._workspace.getNumberHistograms()):
                if spectrum_info.isMonitor(ws_index):
                    monitor_index_list.append(ws_index)
                else:
                    break
            self._n_monitors = len(monitor_index_list)
        return self._n_monitors
