# Unit test for drtsans.mono.spice_data
import pytest
from drtsans.mono.spice_data import SpiceRun


def test_spice_data_constructor():
    """Test constructor and properties access of SpiceRun"""
    # regular constructor
    cg2_run = SpiceRun(beam_line="CG2", ipts_number=828, exp_number=280, scan_number=5, pt_number=1)
    assert cg2_run.beam_line == "CG2"
    assert cg2_run.ipts_number == 828
    assert cg2_run.exp_number == 280
    assert cg2_run.scan_number == 5
    assert cg2_run.pt_number == 1

    # check unique run number generator
    run_number = cg2_run.unique_run_number
    assert run_number == 28000050001


if __name__ == "__main__":
    pytest.main(__file__)
