from __future__ import (absolute_import, division, print_function)
import pytest
from ornl.sans.sns.eqsans.prepare import prepare_direct_beam_center


def test_prepare_direct_beam_center(eqsans_f, eqsans_p):
    x, y = prepare_direct_beam_center(eqsans_f['beamcenter'],
                                      eqsans_p['tubes_to_mask'],
                                      finder_kwargs=dict(CenterX=0.01,
                                                         CenterY=0.02))
    assert (x, y) == pytest.approx((0.0265, 0.0180), abs=1e-04)


if __name__ == '__main__':
    pytest.main()
