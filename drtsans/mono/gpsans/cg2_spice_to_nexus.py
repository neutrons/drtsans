from drtsans.mono.convert_xml_to_nexus import EventNexusConverter


class CG2EventNexusConvert(EventNexusConverter):
    """

    """
    def __init__(self):
        """
        initialization

        work for 48 banks
        """
        super(CG2EventNexusConvert, self).__init__('CG2', 'CG2', 48)

    def get_pid_range(self, bank_id):
        """Set GPSANS bank and pixel ID relation

        Parameters
        ----------
        bank_id: int
            bank ID from 1 to 48

        Returns
        -------
        tuple
            start PID, end PID (assuming PID are consecutive in a bank and end PID is inclusive)

        """
        # Check input valid
        if bank_id < 1 or bank_id > 48:
            raise RuntimeError(f'CG2 (GP-SANS) has 88 banks indexed from 1 to 48. '
                               f'Bank {bank_id} is out of range.')

        # calculate starting PID
        if bank_id <= 24:
            # from 1 to 24: front panel
            start_pid = (bank_id - 1) * 2 * 1024
        else:
            # from 25 to 48: back panel
            start_pid = ((bank_id - 25) * 2 + 1) * 1024

        # calculate end PID
        end_pid = start_pid + 1023

        return start_pid, end_pid
