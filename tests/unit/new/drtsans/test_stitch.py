# import numpy as np
import pytest

r"""
Hyperlinks to drtsans functions
stitch_profiles <https://code.ornl.gov/sns-hfir-scse/sans/sans-backend/blob/next/drtsans/stitch.py>
"""
# from drtsans.stitch import stitch_profiles


def test_stitch(reference_dir):
    r"""
    Test stitching by using a dataset and comparing to expected result.

    Function tested: drtsans.stitch
    Undelying Mantid algorithms:
        Stitch1D https://docs.mantidproject.org/nightly/algorithms/Stitch1D-v3.html

    devs - Jose Borreguero <borreguerojm@ornl.gov>
    SME - Weiren Chen <chenw@ornl.gov>, LiLin He <hel3@ornl.gov>
    """
    pass


if __name__ == '__main__':
    pytest.main([__file__])
