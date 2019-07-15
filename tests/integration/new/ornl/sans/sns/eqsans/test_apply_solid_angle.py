import pytest
from tests.conftest import generate_sans_generic_IDF


@pytest.mark.parametrize('generate_sans_generic_IDF',
        [{'Nx':3, 'Ny': 3, 'dx': 0.00425, 'dy': 0.0055}],
                          indirect=True)
def test_solid_angle(generate_sans_generic_IDF):
    assert generate_sans_generic_IDF == ''


if __name__ == '__main__':
    pytest.main()
