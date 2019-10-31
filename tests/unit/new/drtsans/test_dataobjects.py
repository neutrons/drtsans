from drtsans.dataobjects import IQazimuthal, IQmod
import pytest


def test_IQmod_creation():
    # these are expected to work
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


def test_IQazimuthal_1d_creation():
    # these are expected to work
    IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12])
    IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18])
    IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21])

    # arrays are not parallel
    with pytest.raises(TypeError):
        IQazimuthal([1, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 15], [16, 17, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 18], [19, 20, 21])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 21])

    # not enough arguments
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], wavelength=[16, 17, 18])
    with pytest.raises(TypeError):
        IQazimuthal([1, 2, 3], [4, 5, 6], [7, 8, 9])


def test_IQazimuthal_2d_creation():
    # these are expected to work
    IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]])
    IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                [[13, 14], [15, 16]], [[16, 17], [18, 19]])
    IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])

    # arrays are not parallel
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                    [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                    [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8]], [[10, 11], [12, 13]],
                    [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11]],
                    [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                    [[13, 14]], [[16, 17], [18, 19]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                    [[13, 14], [15, 16]], [[16, 17]], [[19, 20], [21, 22]])
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]], [[10, 11], [12, 13]],
                    [[13, 14], [15, 16]], [[16, 17], [18, 19]], [[19, 20]])

    # not enough arguments
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2], [3, 4]], [[4, 5], [6, 7]], [[7, 8], [9, 10]])

    # qx and qy are linear
    IQazimuthal([[1, 2, 3], [3, 4, 5]], [[4, 5, 6], [6, 7, 8]], [7, 8], [10, 11, 12])
    # qx and qy are linear and not right dimension
    with pytest.raises(TypeError):
        IQazimuthal([[1, 2, 3], [3, 4, 5]], [[4, 5, 6], [6, 7, 8]], [7, 8, 9], [10, 11])


def test_IQcrystal_creation():
    pass


if __name__ == '__main__':
    pytest.main([__file__])
