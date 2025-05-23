from pathlib import Path
import pytest
import os
import subprocess


@pytest.mark.mount_eqsans
@pytest.mark.datarepo
def test_grasp_cg2(reference_dir, has_sns_mount, datarepo_dir):
    if not has_sns_mount:
        pytest.skip("Do not have /SNS properly mounted on this system")

    # arguments
    input_filename = "CG2_9177.nxs.h5"
    input_path = str(Path(datarepo_dir.gpsans) / input_filename)

    output_dir = "test_grasp"
    output_filename = "GPSANS_0009177.dat"
    output_path = os.path.join(os.getcwd(), output_dir, output_filename)

    # run the grasp export command
    subprocess.run(["grasp_cg2", input_path, output_dir], check=False)

    # assert the grasp output file was created
    assert os.path.exists(output_path)

    # used to assert the number of lines in the created file
    line_count = 0

    # open the output grasp file to do some checks
    with open(output_path, "r") as f:
        for line in f:
            # assert first line is correct
            if line_count == 0:
                assert line.startswith("# Panel Scan 13 of 60")
            line_count += 1

    # assert number of lines is correct
    assert line_count == 201

    # cleanup
    os.remove(output_path)
    os.rmdir(output_dir)


if __name__ == "__main__":
    pytest.main([__file__])
