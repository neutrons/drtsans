import pytest

from drtsans.eventslice import resolve_slicing


def test_resolve_slicing():
    reduction_input = {
        "sample": {"runNumber": "12345"},
        "configuration": {"useTimeSlice": False, "useLogSlice": False, "polarization": {"level": "half"}},
    }
    assert resolve_slicing(reduction_input) == (False, False, True)

    with pytest.raises(ValueError) as except_info:
        reduction_input["configuration"] = {
            "useTimeSlice": True,
            "useLogSlice": True,
            "polarization": {"level": "half"},
        }
        resolve_slicing(reduction_input)
    assert "Can't do both time and log slicing" in str(except_info.value)

    reduction_input = {
        "sample": {"runNumber": "1,2,3"},
        "configuration": {"useTimeSlice": True, "useLogSlice": False, "polarization": {"level": "half"}},
    }
    with pytest.raises(ValueError) as except_info:
        resolve_slicing(reduction_input)
    assert "Can't do slicing on summed data sets" in str(except_info.value)

    reduction_input = {
        "sample": {"runNumber": "123"},
        "configuration": {"useTimeSlice": True, "useLogSlice": False, "polarization": {"level": "half"}},
    }
    with pytest.raises(NotImplementedError) as except_info:
        resolve_slicing(reduction_input)
    assert "Time or log slicing on polarized data sets is not implemented yet" in str(except_info.value)

    reduction_input = {
        "configuration": {"useTimeSlice": False, "useLogSlice": True, "polarization": {"level": "half"}},
        "sample": {"runNumber": "123"},
    }
    with pytest.raises(NotImplementedError) as except_info:
        resolve_slicing(reduction_input)
    assert "Time or log slicing on polarized data sets is not implemented yet" in str(except_info.value)
