from drtsans.dataobjects import IQmod
import pytest


def test_IQmod_creation():
    IQmod([1, 2, 3], [4, 5, 6], [7, 8, 9])
    IQmod([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12])
    IQmod([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15])
    IQmod([1, 2, 3], [4, 5, 6], [7, 8, 9], wavelength=[13, 14, 15])

    # intensity isn't 1d
    with pytest.raises(TypeError):
        IQmod([[1, 2], [3, 4]], [4, 5, 6], [7, 8, 9])

    # arrays are not parallel
    with pytest.raises(TypeError):
        IQmod([1, 3], [4, 5, 6], [7, 8, 9])
    with pytest.raises(TypeError):
        IQmod([1, 2, 3], [4, 6], [7, 8, 9])
    with pytest.raises(TypeError):
        IQmod([1, 2, 3], [4, 5, 6], [7, 9])

    # not enough arguments
    with pytest.raises(TypeError):
        IQmod([1, 2, 3], [4, 5, 6])


def test_IQazimuthal():
    pass


def test_IQcrystal():
    pass


if __name__ == '__main__':
    pytest.main([__file__])
