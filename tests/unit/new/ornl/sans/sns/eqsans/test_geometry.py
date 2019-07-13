import pytest
from pytest import approx
from mantid.simpleapi import LoadEventNexus
from ornl.settings import amend_config
from ornl.sans.sns.eqsans.geometry import (sample_aperture_diameter,
                                           source_aperture_diameter)
from ornl.sans.samplelogs import SampleLogs


def test_sample_aperture_diameter(refd):
    with amend_config(data_dir=refd.new.eqsans):
        ws = LoadEventNexus('EQSANS_92353')
    sad = sample_aperture_diameter(ws)
    assert sad == approx(10)
    sad = SampleLogs(ws).single_value('sample-aperture-diameter')
    assert sad == approx(10)


def test_source_aperture_diameter(refd):
    with amend_config(data_dir=refd.new.eqsans):
        ws = LoadEventNexus('EQSANS_92353')
    sad = source_aperture_diameter(ws)
    assert sad == approx(20)
    sad = SampleLogs(ws).single_value('source-aperture-diameter')
    assert sad == approx(20)


if __name__ == '__main__':
    pytest.main()
