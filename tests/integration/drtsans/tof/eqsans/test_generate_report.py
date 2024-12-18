import pytest
import os
import subprocess


@pytest.mark.datarepo
def test_generate_report(datarepo_dir):
    TEST_FILE = os.path.join(datarepo_dir.eqsans, "MC-IL-CO2q_2_reduction_log.hdf")

    cmd = ["generate_report", TEST_FILE]

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    pattern = """drtsan version                                                   1.10.2+d20231019
mantid version                                                              6.8.0
background run number                                                       90335
background_transmission                                        0.9254989091541436
background_transmission error                               0.0029664327599940666
beam center run number                                                      90277
beam center x                                               0.0062761730094059885
beam center y                                               -0.019878669581987165
sample transmission                                            0.8939763002876431
sample transmission error                                    0.002889993726769822
sample run number                                                           90334
"""
    assert result.returncode == 0
    assert result.stderr == ""
    output_list = result.stdout.split("\n")
    pattern_list = pattern.split("\n")
    assert output_list[2:] == pattern_list
