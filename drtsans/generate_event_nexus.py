# Create Event NeXus file
from drtsans.h5_buffer import HDFNode


class BankNode(HDFNode):
    """Node for bank

    """
    def __init__(self):
        """

        """
        super(BankNode, self).__init__()


class EventNeXusWriter(object):
    """
    Write an Event NeXus file
    """
    def __init__(self):
        """ Initialization
        """
        # Bank of events
        self._banks_dict = dict()

        # Meta data
        self._meta_data_dict = dict()

        # Run start time
        self._run_start = None

    def set_counts(self, bank_id, counts, detector_ids):
        self._banks_dict[bank_id] = counts, detector_ids

    def set_meta_data(self, meta_name, value, unit):
        self._meta_data_dict[meta_name] = value, unit

    def set_run_start_time(self, run_start_time):
        self._run_start = run_start_time
