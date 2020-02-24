import os
import pytest


@pytest.mark.parametrize('root',
                         ['/HFIR/CG3', '/SNS/EQSANS', '/HFIR/CG2'],
                         ids=['BIOSANS', 'EQSANS', 'GPSANS'])
def test_mounts(root):
    path = os.path.join(root, 'shared')
    assert os.path.isdir(path)


if __name__ == '__main__':
    pytest.main([__file__])
