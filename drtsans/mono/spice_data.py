# Module containing multiple classes to work with SPICE data
from typing import NamedTuple
import os


class SpiceRun(NamedTuple):
    # Beam line
    beam_line: str
    # IPTS number
    ipts_number: int
    # Experiment number
    exp_number: int
    # Scan number
    scan_number: int
    # Pt number
    pt_number: int

    @property
    def hfir_ipts_dir(self):
        """Standard HFIR SPICE data directory

        Returns
        -------
        str
            Path to IPTS directory

        """
        return f'/HFIR/{self.beam_line}/IPTS-{self.ipts_number}/'

    def locate_spice_file(self, data_dir=None, raise_if_not_exist=True):
        """Locate SPICE file

        Parameters
        ----------
        data_dir : str or None
            None is by default
        raise_if_not_exist: bool
            raise RuntimeError if the file does not exist

        Returns
        -------
        str
            Path to the SPICE file

        """
        # standard spice file name
        spice_file_name = f'{self.beam_line}_exp{self.exp_number}_scan{self.scan_number:04}_{self.pt_number:04}.xml'

        # data file path
        if data_dir is None:
            # default: on the HFIR data server
            data_dir = os.path.join(self.hfir_ipts_dir, f'exp{self.exp_number}/Datafiles')

        spice_file_path = os.path.join(data_dir, spice_file_name)

        # check file existence
        if raise_if_not_exist and not os.path.exists(spice_file_path):
            raise RuntimeError(f'SPICE file {spice_file_name} cannot be found in directory {data_dir}')

        return spice_file_path

    @property
    def unique_run_number(self):
        """Create a unique run number

        Note: HFIR experiment has unique experiment number.

        Returns
        -------
        int
            run number

        """
        return int(f'{self.exp_number}{self.scan_number:04}{self.pt_number:04}')

    def unique_nexus_name(self, nexus_dir=None, raise_if_not_exist=False):
        """

        Parameters
        ----------
        nexus_dir: str, None
            Path to nexus file.  If None
        raise_if_not_exist: bool
            Check the nexus file exists or not.  Raise exception if it does not

        Returns
        -------
        str
            full path to Nexus file

        """
        # standard base name
        base_nexus_name = f'{self.beam_line}_{self.unique_run_number:012}.nxs.h5'

        if nexus_dir is None:
            # use standard defalt
            nexus_dir = os.path.join(self.hfir_ipts_dir, f'shared/Exp{self.exp_number}')

        nexus_path = os.path.join(nexus_dir, base_nexus_name)

        if raise_if_not_exist and not os.path.exists(nexus_path):
            raise RuntimeError(f'Spice converted Nexus file {base_nexus_name} does not exist in '
                               f'directory {nexus_dir}')

        return nexus_path
